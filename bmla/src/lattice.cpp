#include "ixn_param.hpp"
#include <iostream>
#include <fstream>
#include <numeric>
#include "../include/general.hpp"
#include "math.h"
#include <ctime>
#include <sstream>
#include <random>

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Struct to hold a lattice Site
	****************************************/

	// Constructor
	Site::Site() : Site(0,0,0,nullptr) { dim=0; };
	Site::Site(int xIn) : Site(xIn,0,0,nullptr) { dim=1; };
	Site::Site(int xIn, Species *spIn) : Site(xIn,0,0,spIn) { dim=1; };
	Site::Site(int xIn, int yIn) : Site(xIn,yIn,0,nullptr) { dim=2; };
	Site::Site(int xIn, int yIn, Species *spIn) : Site(xIn,yIn,0,spIn) { dim=2; };
	Site::Site(int xIn, int yIn, int zIn) : Site(xIn, yIn, zIn, nullptr) { dim=3; };
	Site::Site(int xIn, int yIn, int zIn, Species *spIn) {
		dim=3;
		x = xIn;
		y = yIn;
		z = zIn;
		sp = spIn;
		binary = true;
		_prob_empty = 0.0;
	};	
	Site::Site(const Site& other) {
		_copy(other);
	};
	Site::Site(Site&& other) {
		_copy(other);
		other._reset();
	};
	Site& Site::operator=(const Site& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	Site& Site::operator=(Site&& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
	};
	Site::~Site() {
		_clean_up();
	};
	void Site::_clean_up() {
		// Nothing....
	};
	void Site::_reset() {
		dim = 0;
		x = 0;
		y = 0;
		z = 0;
		sp = nullptr;
		nbrs.clear();
		hidden_conns.clear();
		binary = true;
		_prob_empty = 0.0;
		_probs.clear();
	};
	void Site::_copy(const Site& other) {
		dim = other.dim;
		x = other.x;
		y = other.y;
		z = other.z;
		sp = other.sp;
		nbrs = other.nbrs;
		hidden_conns = other.hidden_conns;
		binary = other.binary;
		_prob_empty = other._prob_empty;
		_probs = other._probs;
	};


	// Comparator
	bool operator <(const Site& a, const Site& b) {
		if (a.dim == 1) {
	    	return a.x < b.x;
	    } else if (a.dim == 2) {
	    	return std::tie(a.x, a.y) < std::tie(b.x, b.y);
	    } else if (a.dim == 3) {
	    	return std::tie(a.x, a.y, a.z) < std::tie(b.x, b.y, b.z);
	    } else {
	    	return false;
	    };
	};
	bool operator==(const Site& a, const Site& b) {
		if (a.dim == 1) {
	    	return a.x == b.x;
	    } else if (a.dim == 2) {
	    	return std::tie(a.x, a.y) == std::tie(b.x, b.y);
	    } else if (a.dim == 3) {
	    	return std::tie(a.x, a.y, a.z) == std::tie(b.x, b.y, b.z);
	    } else {
	    	return false;
	    };
	}; 
	std::ostream& operator<<(std::ostream& os, const Site& s)
	{
		if (s.dim == 1) {
			if (s.sp) {
			    return os << s.x << " " << s.sp->name();
			} else {
			    return os << s.x << " empty";
			};	    
		} else if (s.dim == 2) {
			if (s.sp) {
			    return os << s.x << " " << s.y << " " << s.sp->name();
			} else {
			    return os << s.x << " " << s.y << " empty";
			};	    
		} else if (s.dim == 3) {
			if (s.sp) {
			    return os << s.x << " " << s.y << " " << s.z << " " << s.sp->name();
			} else {
			    return os << s.x << " " << s.y << " " << s.z << " empty";
			};
	    } else {
	    	return os;
	    };
	};

	void Site::normalize_probs() {
		double tot = _prob_empty;
		for (auto mp: _probs) {
			tot += mp.second;
		};
		_prob_empty /= tot;
		for (auto mp: _probs) {
			_probs[mp.first] /= tot;
		};
	};

	// Get a probability
	// nullptr for empty
	double Site::get_prob(Species *sp) {
		// nullptr for prob of empty
		if (sp == nullptr) {
			return _prob_empty;
		};

		auto it = _probs.find(sp);
		if (it != _probs.end()) {
			return it->second;
		} else {
			return 0.0;
		};
	};
	void Site::set_prob(Species *sp, double val) {
		// nullptr for prob of empty
		if (sp == nullptr) {
			_prob_empty = val;
			return;
		};

		// Store
		_probs[sp] = val;

		// Increment counts on species
		// Counts
		sp->count_increment(val);
		// NNs
		double p;
		// Go through neighbor sites
		for (auto nbr_it: nbrs) {
			nbr_it->mult_probs_as_nbr(sp,val);
		};
	};
	void Site::mult_probs_as_nbr(Species *sp, double val) {
		for (auto mp: _probs) {
			mp.first->nn_count_increment(sp,val*mp.second);
			if (mp.first != sp) {
				sp->nn_count_increment(mp.first, val*mp.second);
			};
		};
	};


	/****************************************
	Lattice
	****************************************/

	/********************
	Constructor
	********************/

	Lattice::Lattice(int dim, int box_length)
	{
		if (dim != 1 && dim != 2 && dim != 3) {
			std::cerr << "ERROR: only dimensions 1,2,3 are supported for Lattice." << std::endl;
			exit(EXIT_FAILURE);
		};
		_dim = dim;
		_box_length = box_length;
		_hidden_layer_exists = false;

		// Make a fully linked list of sites
		if (dim == 1) {
			for (int x=1; x<=box_length; x++) {
				_latt.push_back(Site(x));					
			};
		} else if (dim == 2) {
			for (int x=1; x<=box_length; x++) {
				for (int y=1; y<=box_length; y++) {
					_latt.push_back(Site(x,y));					
				};
			};
		} else if (dim == 3) {
			for (int x=1; x<=box_length; x++) {
				for (int y=1; y<=box_length; y++) {
					for (int z=1; z<=box_length; z++) {
						_latt.push_back(Site(x,y,z));					
					};
				};
			};
		};

		// Set up neighbors
		std::vector<Site> nbrs;
		latt_it lit = _latt.begin();
		while (lit != _latt.end()) {
			// Neighbors
			nbrs.clear();
			if (lit->x-1 >= 1) {
				nbrs.push_back(Site(lit->x-1,lit->y,lit->z));
			};
			if (lit->x+1 <= box_length) {
				nbrs.push_back(Site(lit->x+1,lit->y,lit->z));
			};
			if (dim == 2 || dim == 3) {
				if (lit->y-1 >= 1) {
					nbrs.push_back(Site(lit->x,lit->y-1,lit->z));
				};
				if (lit->y+1 <= box_length) {
					nbrs.push_back(Site(lit->x,lit->y+1,lit->z));
				};
			};
			if (dim == 3) {
				if (lit->z-1 >= 1) {
					nbrs.push_back(Site(lit->x,lit->y,lit->z-1));
				};
				if (lit->z+1 <= box_length) {
					nbrs.push_back(Site(lit->x,lit->y,lit->z+1));
				};
			};

			// Go through neighbors
			for (auto nbr: nbrs) {
				// Add as nbr
				if (dim == 1) {
					lit->nbrs.push_back(_look_up(nbr.x));
				} else if (dim == 2) {
					lit->nbrs.push_back(_look_up(nbr.x,nbr.y));
				} else if (dim == 3) {
					lit->nbrs.push_back(_look_up(nbr.x,nbr.y,nbr.z));
				};
			};

			// Next
			lit++;
		};

	};
	Lattice::Lattice() {
		_dim = 0;
		_box_length = 0;
	};
	Lattice::Lattice(const Lattice& other) {
		_copy(other);
	};
	Lattice::Lattice(Lattice&& other) {
		_copy(other);
		other._reset();
	};
	Lattice& Lattice::operator=(const Lattice& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	Lattice& Lattice::operator=(Lattice&& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
	};
	Lattice::~Lattice() {
		_clean_up();
	};

	void Lattice::_clean_up() {
		// Nothing...
	};
	void Lattice::_copy(const Lattice& other) {
		_latt = other._latt;
		_box_length = other._box_length;
		_sp_map = other._sp_map;
		_sp_vec = other._sp_vec;
		_hidden_layer_exists = other._hidden_layer_exists;
	};
	void Lattice::_reset() {
		_box_length = 0;
		_sp_map.clear();
		_sp_vec.clear();
		_latt.clear();
		_hidden_layer_exists = false;
	};

	/********************
	Getters
	********************/

	int Lattice::dim() const {
		return _dim;
	};

	/********************
	Add a species
	********************/

	void Lattice::add_species(Species *sp) {
		if (sp) {
			_sp_map[sp->name()] = sp;
			_sp_vec.push_back(sp);
		};
	};

	/********************
	Indicate that the hidden unit exists
	********************/

	void Lattice::set_hidden_layer_exists() {
		_hidden_layer_exists = true;
	};

	/********************
	Find a pointer to a site by index
	********************/

	Site* Lattice::get_site(int x) {
		if (_dim != 1) {
			std::cerr << "ERROR: dim wrong in get_site" << std::endl;
			exit(EXIT_FAILURE);
		};
		return &(*(_look_up(x)));
	};
	Site* Lattice::get_site(int x, int y) {
		if (_dim != 2) {
			std::cerr << "ERROR: dim wrong in get_site" << std::endl;
			exit(EXIT_FAILURE);
		};
		return &(*(_look_up(x,y)));
	};
	Site* Lattice::get_site(int x, int y, int z) {
		if (_dim != 3) {
			std::cerr << "ERROR: dim wrong in get_site" << std::endl;
			exit(EXIT_FAILURE);
		};
		return &(*(_look_up(x,y,z)));
	};

	/********************
	Clear, size
	********************/

	void Lattice::clear() { 
		// Reset the counts and nns
		for (auto sp: _sp_vec) {
			sp->reset_counts();
		};

		// Clear the lattice
		for (auto it = _latt.begin(); it != _latt.end(); it++) {
			it->sp = nullptr;
		};
	};
	int Lattice::size() { return _latt.size(); };

	/********************
	Make a mol
	********************/

	bool Lattice::make_mol(latt_it s, Species *sp) {
		// Check not already occupide
		if (s->sp) {
			std::cerr << "ERROR site already occupied" << std::endl;
			return false;
		};

		// Make mol at empty site
		return make_mol_at_empty(s,sp);
	};

	bool Lattice::replace_mol(latt_it s, Species *sp) {
		// Check not the same as already there
		if (s->sp == sp) {
			return true;
		};

		// If already occupied, remove
		if (s->sp) {
			// Erase mol
			erase_mol(s);
		};

		// Make mol at empty site
		return make_mol_at_empty(s,sp);
	};

	bool Lattice::make_mol_at_empty(latt_it s, Species *sp) {
		/*
		std::cout << "Making: " << s->x << " " << s->y << " " << s->z << std::flush;
		int ctr=0;
		for (auto nbr_it: s->nbrs) { // go through nbrs
			if (nbr_it->sp) { // is nbr site occ
				ctr++;
			};
		};
		std::cout << " no nbrs: " << ctr << std::flush;
		std::cout << " count now: " << sp->nn_count(_sp_map["A"]) << std::flush;
		*/

		// Make
		s->sp = sp;

		// Update count and nns
		increment_species_counts_binary(s,sp,1.0);

		// std::cout << " now: count: " << sp->nn_count[_sp_map["A"]] << std::endl;

		return true;
	};

	// Increment counts on species
	void Lattice::increment_species_counts_binary(latt_it s, Species *sp, double inc) {
		sp->count_increment(inc);
		for (auto nbr_it: s->nbrs) { // go through nbrs
			if (nbr_it->sp) { // is nbr site occ
				// Update bidirectionally, unless it's the same species
				sp->nn_count_increment(nbr_it->sp,inc);
				if (nbr_it->sp != sp) {
					nbr_it->sp->nn_count_increment(sp,inc);
				};
			};
		};
	};

	/********************
	Erase a mol
	********************/

	bool Lattice::erase_mol(latt_it s)
	{
		// Check not empty
		if (!(s->sp)) {
			std::cerr << "ERROR site not occupied" << std::endl;
			return false;
		};

		/*
		std::cout << "Erasing: " << s->x << " " << s->y << " " << s->z;
		int ctr=0;
		for (auto nbr_it: s->nbrs) { // go through nbrs
			if (nbr_it->sp) { // is nbr site occ
				ctr++;
			};
		};
		std::cout << " no nbrs: " << ctr << " count now: " << s->sp->nn_count[_sp_map["A"]];
		*/

		// Update count and nns
		increment_species_counts_binary(s,s->sp,-1.0);

		// std::cout << " now: count: " << s->sp->nn_count[_sp_map["A"]] << std::endl;

		// Erase
		s->sp = nullptr;

		return true;
	};

	/********************
	Write lattice to a file
	********************/

	void Lattice::write_to_file(std::string fname) 
	{
		std::ofstream f;
		f.open (fname);
		for (auto l: _latt) {
			if (l.sp) {
				f << l << "\n";
			};
		};
		f.close();
	};

	/********************
	Read lattice from a file
	********************/

	void Lattice::read_from_file(std::string fname)
	{
		// Clear the current lattice
		clear();

		std::ifstream f;
		f.open(fname);
		std::string x="",y="",z="";
		std::string sp="";
		std::string line;
		std::istringstream iss;
		if (f.is_open()) { // make sure we found it
			while (getline(f,line)) {
				if (line == "") { continue; };
				iss = std::istringstream(line);
			    iss >> x;
			    if (_dim == 2 || _dim == 3) {
				    iss >> y;
			    };
			    if (_dim == 3) {
				    iss >> z;
			    };
			    iss >> sp;
		    	// Add to map
		    	if (_dim == 1) {
			    	make_mol(_look_up(atoi(x.c_str())),_sp_map[sp]);
			    } else if (_dim == 2) {
			    	make_mol(_look_up(atoi(x.c_str()),atoi(y.c_str())),_sp_map[sp]);
			    } else if (_dim == 3) {
			    	make_mol(_look_up(atoi(x.c_str()),atoi(y.c_str()),atoi(z.c_str())),_sp_map[sp]);
			    };
		    	sp=""; x=""; y=""; z="";
			};
		};
		f.close();

		// std::cout << "Read: " << fname << std::endl;
	};

	/********************
	Sample
	********************/

	void Lattice::sample(bool binary) {

		// Declarations
		latt_it it;
		double energy;
		std::map<Species*, std::vector<HiddenUnit*>>::iterator it_hups;
		std::vector<double> probs; // probabilities
		std::vector<double> props; // propensities
		int i_chosen;
		latt_it lit;

		// Copy lattice
		lattice latt_new(_latt);

		// Clear everything existing
		clear();

		// Go through all lattice sites
		for (it = latt_new.begin(); it != latt_new.end(); it++) {

			// Look up the site in our latt
			if (_dim==1) {
				lit = _look_up(it->x);
			} else if (_dim == 2) {
				lit = _look_up(it->x,it->y);
			} else if (_dim == 3) {
				lit = _look_up(it->x,it->y,it->z);
			};

			// Clear propensities
			props.clear();
			props.push_back(0.0);

			// Propensity for no spin is exp(0) = 1
			props.push_back(1.0);

			// Prob
			probs.clear();
			probs.push_back(1.0);

			// Go through all possible species this could be, calculate propensities
			for (auto sp_new: _sp_vec) {
				// Bias
				energy = sp_new->h();

				// NNs for J
				for (auto it_nbr: it->nbrs) {
					// Occupied?
					if (it_nbr->sp) {
						energy += sp_new->j(it_nbr->sp);
					};
				};

				// Hidden layer exists?
				if (_hidden_layer_exists) {

					// Check if this species has connections to hidden units
					it_hups = it->hidden_conns.find(sp_new);
					if (it_hups != it->hidden_conns.end()) {
						// Yes it does - sum them up! Go over hidden units
						for (auto hup: it_hups->second) {
							// Weight * value of spin (0 or 1 if binary, else a prob)
							energy += sp_new->w() * hup->get();
						};
					};
				};

				// Append prob
				energy = exp(energy);
				probs.push_back(energy);
				props.push_back(props.back()+energy);

			};

			// Binarize if needed
			if (binary) {

				// Binary
				lit->binary = true;

				// Sample RV
				// i_chosen = sample_prob_vec(probs);
				i_chosen = sample_prop_vec(props);

				if (i_chosen==0) {
					// Flip down (new spin = 0)
					// Already guaranteed empty
				} else {
					// Make the appropriate species at this site (guaranteed empty)
					make_mol_at_empty(lit,_sp_vec[i_chosen-1]);
				};
			} else {

				// Probabilistic
				lit->binary = false;

				// Normalize the probs
				double tot=0.0;
				for (auto t: probs) {
					tot += t;
				};

				// Write into species
				lit->set_prob(nullptr,probs[0]/tot);
				for (int i=1; i<probs.size(); i++) {
					lit->set_prob(_sp_vec[i-1],probs[i]/tot);
				};
				/*
				for (auto pr: probs) {
					std::cout << pr/tot << " ";
				};
				std::cout << std::endl;
				*/
			};
		};
	};

	/********************
	Sample probabilities/propensities
	********************/

	// Sample an unnormalized probability vector
	int Lattice::sample_prob_vec(std::vector<double> &probs) {
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator(seed);
		std::discrete_distribution<int> dd(probs.begin(),probs.end());
		return dd(generator);
	};

	// Sample a vector of propensities (cumulative probabilities)
	int Lattice::sample_prop_vec(std::vector<double> &props) {
		// Sample RV
		double r = randD(0.0,props.back());

		// Find interval
		for (int i=0; i<props.size()-1; i++) {
			if (props[i] <= r && r <= props[i+1]) {
				return i;
			};
		};
		return 0; // never get here
	};

	/****************************************
	Lattice - PRIVATE
	****************************************/

	/********************
	Lookup a site iterator from x,y,z
	********************/

	latt_it Lattice::_look_up(int x) {
		// Figure out index in list
		int n = x-1;

		// Grab
		latt_it it = _latt.begin();
		std::advance(it,n);
		return it;
	};
	latt_it Lattice::_look_up(int x, int y) {
		// Figure out index in list
		int n = (x-1)*_box_length + y-1;

		// Grab
		latt_it it = _latt.begin();
		std::advance(it,n);
		return it;
	};
	latt_it Lattice::_look_up(int x, int y, int z) {
		// Figure out index in list
		int n = (x-1)*_box_length*_box_length + (y-1)*_box_length + z-1;

		// Grab
		latt_it it = _latt.begin();
		std::advance(it,n);
		return it;
	};




};