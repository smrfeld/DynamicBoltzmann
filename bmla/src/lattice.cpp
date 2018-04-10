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
		_sp_binary = spIn;
		_binary = true; // default
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
		_sp_binary = nullptr;
		nbrs.clear();
		hidden_conns.clear();
		_binary = true;
		_prob_empty = 0.0;
		_probs.clear();
	};
	void Site::_copy(const Site& other) {
		dim = other.dim;
		x = other.x;
		y = other.y;
		z = other.z;
		_sp_binary = other._sp_binary;
		nbrs = other.nbrs;
		hidden_conns = other.hidden_conns;
		_binary = other._binary;
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
		Species *sp;
		if (s.binary()) {
			sp = s.get_species_binary();
			if (sp)	{
				if (s.dim == 1) {
				    return os << s.x << " " << sp->name();
				} else if (s.dim == 2) {
				    return os << s.x << " " << s.y << " " << sp->name();
				} else if (s.dim == 3) {
				    return os << s.x << " " << s.y << " " << s.z << " " << sp->name();
				};
			};
		};
    	
    	return os;
	};

	// Get a probability
	// nullptr for empty
	double Site::get_prob(Species *sp) const {
		// Binary?
		if (_binary) {
			// Binary
			if (sp == _sp_binary) {
				return 1.0;
			} else {
				return 0.0;
			};
		} else {
			// Probabilistic
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
	};

	// Get all probs - if binary, returns an empty
	const std::map<Species*, double>& Site::get_probs() const {
		return _probs;
	};

	// Get J and K activations
	// Go through all possible species probabilities x J of the coupling for the given species
	double Site::get_act_j(Species *sp) const {
		double act=0.0;
		// Go through all nbrs
		for (auto lit: nbrs) {
			// Get all probs
			const std::map<Species*, double> prs = lit->get_probs();
			// Go through all probs
			for (auto pr: prs) {
				// J * prob
				act += sp->j(pr.first) * pr.second;
			};
		};
		return act;
	};
	double Site::get_act_k(Species *sp) const {
		return 0.0;
	};

	// Set probability
	// Pass nullptr to set probability of being empty
	void Site::set_prob(Species *sp, double prob) {
		// Check if currently binary
		if (_binary) {
			_clear_binary();
		};

		// nullptr for prob of empty
		if (sp == nullptr) {
			_prob_empty = prob;
			return;
		};

		// Store
		_probs[sp] = prob;

		// Increment counts on species
		_add_counts_on_species_prob(sp,prob);
	};

	// Increment NN counts given that we are neighboring this species with this prob
	void Site::increment_nn_counts_for_neighbor_prob(Species *sp_nbr, double prob) {
		// Go through all species/probs
		for (auto mp: _probs) {
			// Increment counts bidirectionally
			mp.first->nn_count_increment(sp_nbr,prob*mp.second);
			if (mp.first != sp_nbr) {
				sp_nbr->nn_count_increment(mp.first, prob*mp.second);
			};
		};
	};

	// If binary, gets the species
	// Else returns nullptr if empty OR probabilistic
	Species* Site::get_species_binary() const {
		if (_binary) {
			return _sp_binary; 
		} else {
			return nullptr;
		};
	};

	// Sets a site to binary with some species
	// Pass nullptr to make an empty binary site
	void Site::set_species_binary(Species *sp) {
		// Check if currently probabilistic
		if (!_binary) {
			_clear_prob();
		};

		// Is already occupied? Then remove counts
		if (_sp_binary) {
			_remove_counts_on_species_binary(_sp_binary);
		};

		// Set
		_sp_binary = sp;
		
		// Increment counts
		_add_counts_on_species_binary(_sp_binary);
	};

	// Set a site to be empty
	void Site::set_site_empty_binary() {
		// Check if currently probabilistic
		if (!_binary) {
			_clear_prob();
		};

		// Check if already empty
		if (!_sp_binary) {
			return;
		};

		// Remove counts on existing species
		_remove_counts_on_species_binary(_sp_binary);

		// Clear species
		_sp_binary = nullptr;
	};

	// Check if the site is binary
	bool Site::binary() const { return _binary; };

	// Add/Remove counts on a given species (not _sp_binary, unless it is passed)
	void Site::_remove_counts_on_species_binary(Species *sp) {
		// Count
		sp->count_increment(-1.0);
		// NN
		Species *sp_nbr = nullptr;
		for (auto nbr_it: nbrs) { // go through nbrs
			sp_nbr = nbr_it->get_species_binary();
			if (sp_nbr) { // is nbr site occ
				// Update bidirectionally, unless it's the same species
				sp->nn_count_increment(sp_nbr,-1.0);
				if (sp_nbr != sp) {
					sp_nbr->nn_count_increment(sp,-1.0);
				};
			};
		};
	};
	void Site::_add_counts_on_species_binary(Species *sp) {
		// Count
		sp->count_increment(1.0);
		// NN
		Species *sp_nbr = nullptr;
		for (auto nbr_it: nbrs) { // go through nbrs
			sp_nbr = nbr_it->get_species_binary();
			if (sp_nbr) { // is nbr site occ
				// Update bidirectionally, unless it's the same species
				sp->nn_count_increment(sp_nbr,1.0);
				if (sp_nbr != sp) {
					sp_nbr->nn_count_increment(sp,1.0);
				};
			};
		};
	};
	void Site::_remove_counts_on_species_prob(Species *sp, double prob) {
		// Counts
		sp->count_increment(-1.0*prob);
		// NNs
		for (auto nbr_it: nbrs) {
			nbr_it->increment_nn_counts_for_neighbor_prob(sp,-1.0*prob);
		};
	};
	void Site::_add_counts_on_species_prob(Species *sp, double prob) {
		// Counts
		sp->count_increment(prob);
		// NNs
		for (auto nbr_it: nbrs) {
			nbr_it->increment_nn_counts_for_neighbor_prob(sp,prob);
		};
	};

	// Clear a site from being binary/probabilistic
	// Note: flips to the opposite mode
	void Site::_clear_binary() {
		// Remove counts if needed
		if (_sp_binary) {
			_remove_counts_on_species_binary(_sp_binary);
		};

		// Remove
		_sp_binary = nullptr;

		// Not binary now
		_binary = false;
	};
	void Site::_clear_prob() {
		// Remove counts if needed
		for (auto pr: _probs) {
			_remove_counts_on_species_prob(pr.first,pr.second);
		};

		// Remove
		_probs.clear();
		_prob_empty = 0.0;

		// Now binary
		_binary = true;

		// Set default binary state = empty
		_sp_binary = nullptr;
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
		// Clear the lattice
		for (auto it = _latt.begin(); it != _latt.end(); it++) {
			it->set_site_empty_binary();
		};

		// Reset the counts and nns
		for (auto sp: _sp_vec) {
			sp->reset_counts();
		};
	};
	int Lattice::size() { return _latt.size(); };

	/********************
	Write lattice to a file
	********************/

	void Lattice::write_to_file(std::string fname) 
	{
		std::ofstream f;
		f.open (fname);
		for (auto l: _latt) {
			f << l << "\n";
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
		latt_it lit;
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
		    	// Add to lattice
		    	if (_dim == 1) {
		    		lit = _look_up(atoi(x.c_str()));
			    } else if (_dim == 2) {
		    		lit = _look_up(atoi(x.c_str()),atoi(y.c_str()));
			    } else if (_dim == 3) {
			    	lit = _look_up(atoi(x.c_str()),atoi(y.c_str()),atoi(z.c_str()));
			    };
	    		lit->set_species_binary(_sp_map[sp]);
	    		// Reset
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
				energy += it->get_act_j(sp_new);

				// Triplets for K
				if (_dim == 1) { // Only dim 1 currently supported
					// ...
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

				// Sample RV
				// i_chosen = _sample_prob_vec(probs);
				i_chosen = _sample_prop_vec(props);

				if (i_chosen==0) {
					// Flip down (new spin = 0)
					// Already guaranteed empty, since lattice was cleared
				} else {
					// Make the appropriate species at this site (guaranteed empty)
					lit->set_species_binary(_sp_vec[i_chosen-1]);
				};
			} else {

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

	/********************
	Sample probabilities/propensities
	********************/

	// Sample an unnormalized probability vector
	int Lattice::_sample_prob_vec(std::vector<double> &probs) {
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator(seed);
		std::discrete_distribution<int> dd(probs.begin(),probs.end());
		return dd(generator);
	};

	// Sample a vector of propensities (cumulative probabilities)
	int Lattice::_sample_prop_vec(std::vector<double> &props) {
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

};