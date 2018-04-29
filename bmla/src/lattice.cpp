#include "ixn_param.hpp" // also includes lattice header
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
	Class to hold a connection from visible to hidden
	****************************************/

	// Constructor
	ConnectionVH::ConnectionVH(Site *site, HiddenUnit *hidden_unit, std::vector<IxnParam*> ips) {
		_site = site;
		_hidden_unit = hidden_unit;
		for (auto ip: ips) {
			add_ixn_param(ip);
		};	
	};
	ConnectionVH::ConnectionVH(const ConnectionVH& other) {
		_copy(other);
	};
	ConnectionVH::ConnectionVH(ConnectionVH&& other) {
		_copy(other);
		other._reset();
	};
	ConnectionVH& ConnectionVH::operator=(const ConnectionVH& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	ConnectionVH& ConnectionVH::operator=(ConnectionVH&& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
	};
	ConnectionVH::~ConnectionVH() {
		_clean_up();
	};
	void ConnectionVH::_clean_up() {
		// Nothing....
	};
	void ConnectionVH::_reset() {
		_site = nullptr;
		_hidden_unit = nullptr;
		_ips_visible_hidden.clear();
		_ips_hidden_visible.clear();
	};
	void ConnectionVH::_copy(const ConnectionVH& other) {
		_site = other._site;
		_hidden_unit = other._hidden_unit;
		_ips_visible_hidden = other._ips_visible_hidden;
		_ips_hidden_visible = other._ips_hidden_visible;
	};

	// Add ixn param
	void ConnectionVH::add_ixn_param(IxnParam* ip) {
		// Get the species associated with this ixn param
		std::vector<SpeciesVH> sp_vec = ip->get_species_conn();

		// Go through the species
		for (auto spvh: sp_vec) {
			// Store both ways
			_ips_visible_hidden[spvh.sp_visible][spvh.sp_hidden].push_back(ip);
			_ips_hidden_visible[spvh.sp_hidden][spvh.sp_visible].push_back(ip);
		};
	};

	// Get for a species on the visible unit
	double ConnectionVH::get_act_visible(Species* sp_visible) {
		double act=0.0;
		auto it = _ips_visible_hidden.find(sp_visible);
		if (it != _ips_visible_hidden.end()) {
			// Go through all hidden species
			for (auto iph: it->second) {
				// Go through all ixn params
				for (auto ip: iph.second) {
					// Weight (from ixn param) * hidden units value for this hidden species
					act += ip->get() * _hidden_unit->get_prob(iph.first);
				};
			};
		};
		return act;
	};

	// Get activation for a hidden
	double ConnectionVH::get_act_hidden(HiddenSpecies* sp_hidden) {
		double act=0.0;
		auto it = _ips_hidden_visible.find(sp_hidden);
		if (it != _ips_hidden_visible.end()) {
			// Go through all visible species
			for (auto ipv: it->second) {
				// Go through all ixn params
				for (auto ip: ipv.second) {
					// Weight (from ixn param) * sites value for this visible species
					act += ip->get() * _site->get_prob(ipv.first);
				};
			};
		};
		return act;
	};

	/****************************************
	Struct to hold a lattice Site
	****************************************/

	/********************
	Constructor
	********************/

	Site::Site(int xIn) : Site(xIn,0,0) { dim=1; };
	Site::Site(int xIn, int yIn) : Site(xIn,yIn,0) { dim=2; };
	Site::Site(int xIn, int yIn, int zIn) {
		dim=3;
		x = xIn;
		y = yIn;
		z = zIn;
		_prob_empty = 1.0; // default = empty
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
		nbrs.clear();
		nbrs_triplets.clear();
		nbrs_quartics.clear();
		_hidden_conns.clear();
		_prob_empty = 0.0;
		_probs.clear();
	};
	void Site::_copy(const Site& other) {
		dim = other.dim;
		x = other.x;
		y = other.y;
		z = other.z;
		nbrs = other.nbrs;
		nbrs_triplets = other.nbrs_triplets;
		nbrs_quartics = other.nbrs_quartics;
		_hidden_conns = other._hidden_conns;
		_prob_empty = other._prob_empty;
		_probs = other._probs;
	};

	/********************
	Comparator
	********************/

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
		if (s.empty()) {
			return os;
		};

		if (s.dim == 1) {
		    os << s.x << " ";
		} else if (s.dim == 2) {
		    os << s.x << " " << s.y << " ";
		} else if (s.dim == 3) {
		    os << s.x << " " << s.y << " " << s.z << " ";
		};

		// Get probs
		const std::map<Species*, double> prs = s.get_probs();

		// All probs
		for (auto pr: prs) {
			if (pr.second > 0.0) {
				os << pr.first->name() << " " << pr.second << " ";	
			};
		};

	    return os;
    };

	/********************
	Add a hidden conn
	********************/

	void Site::add_visible_hidden_conn(ConnectionVH* connvh) {
		_hidden_conns.push_back(connvh);
	};

	/********************
	Add a species possibility
	********************/

	void Site::add_species_possibility(Species* sp) {
		_probs[sp] = 0.0;
	};

	/********************
	Get probability
	********************/

	// nullptr for empty
	double Site::get_prob(Species *sp) const {
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

	// Get all probs (excluding prob of empty)
	const std::map<Species*, double>& Site::get_probs() const {
		return _probs;
	};

	// Set probability
	// Pass nullptr to set probability of being empty
	void Site::set_prob(Species *sp, double prob) {

		// nullptr for prob of empty
		if (sp == nullptr) {
			_prob_empty = prob;
			return;
		};

		// Remove the old counts
		if (_probs[sp] > 0.0) {
			_remove_counts_on_species(sp,_probs[sp]);
		};

		// Store
		// std::cout << "Set prob for site: " << x << " species " << sp->name() << " to: " << prob << std::endl;
		_probs[sp] = prob;

		// Increment counts on species
		_add_counts_on_species(sp,prob);
	};

	// Set a site to be empty
	void Site::set_site_empty() {
		// Remove counts on existing species
		for (auto pr: _probs) {
			if (pr.second > 0.0) {
				_remove_counts_on_species(pr.first,pr.second);
				// Clear
				_probs[pr.first] = 0.0;
			};
		};

		// Empty prob = 1
		_prob_empty = 1.0;
	};

	// Set site to have binary probs
	void Site::set_site_binary(Species *sp) {
		// Clear
		set_site_empty();
		// Make
		set_prob(sp,1.0);
	};

	/********************
	Get activations for a given species
	********************/

	double Site::get_activation(Species *sp) const {
		double energy;

		// Bias
		energy = sp->h();

		// NNs - go through neihbors
		for (auto lit: nbrs) {
			// Get all probs
			const std::map<Species*, double> prs = lit->get_probs();
			// Go through all probs
			for (auto pr: prs) {
				// J * prob
				energy += sp->j(pr.first) * pr.second;
			};
		};

		// Triplets - go through all pairs to consider
		for (auto trip: nbrs_triplets) {
			// Get all probs
			const std::map<Species*, double> prs1 = trip.lit1->get_probs();
			const std::map<Species*, double> prs2 = trip.lit2->get_probs();
			// Go through all probs
			for (auto pr1: prs1) {
				for (auto pr2: prs2) {
					// K * prob * prob
					energy += sp->k(pr1.first,pr2.first) * pr1.second * pr2.second;
				};
			};
		};

		// Conn to hidden layer - go through conns
		for (auto connvh: _hidden_conns) {
			energy += connvh->get_act_visible(sp);
		};

		return energy;
	};

	/********************
	Check if site is empty
	********************/

	bool Site::empty() const {
		if (_prob_empty == 1.0) {
			return true;
		} else {
			return false;
		};
	};

	/********************
	Binarize the set
	********************/

	void Site::binarize() {
		// Propensity vector
		std::vector<double> props;
		std::vector<Species*> sp_vec;
		props.push_back(0.0);
		props.push_back(_prob_empty);
		for (auto pr: _probs) {
			props.push_back(props.back()+pr.second);
			sp_vec.push_back(pr.first);
		};

		int i = sample_prop_vec(props);
		if (i==0) {
			set_site_empty();
		} else {
			set_site_binary(sp_vec[i-1]);
		};
	};

	/****************************************
	Site - PRIVATE
	****************************************/

	void Site::_remove_counts_on_species(Species *sp, double prob) {
		_add_counts_on_species(sp,-1.0*prob);
	};
	void Site::_add_counts_on_species(Species *sp, double prob) {

		// Counts
		sp->count_increment(prob);

		// NNs, if needed
		for (auto nbr_it: nbrs) {
			// Get all the probs
			const std::map<Species*, double> prs = nbr_it->get_probs();
			for (auto pr: prs) {
				// Increment
				sp->nn_count_increment(pr.first,prob*pr.second);
			};
		};

		// Triplets, if needed
		for (auto trip: nbrs_triplets) {
			// Get all the probs
			const std::map<Species*, double> prs1 = trip.lit1->get_probs();
			const std::map<Species*, double> prs2 = trip.lit2->get_probs();
			for (auto pr1: prs1) {
				for (auto pr2: prs2) {
					// Increment
					sp->triplet_count_increment(pr1.first,pr2.first,prob*pr1.second*pr2.second);
				};
			};
		};

		// Quartics, if needed
		for (auto quart: nbrs_quartics) {
			// Get all the probs
			const std::map<Species*, double> prs1 = quart.lit1->get_probs();
			const std::map<Species*, double> prs2 = quart.lit2->get_probs();
			const std::map<Species*, double> prs3 = quart.lit3->get_probs();
			for (auto pr1: prs1) {
				for (auto pr2: prs2) {
					for (auto pr3: prs3) {
						// Increment
						sp->quartic_count_increment(pr1.first,pr2.first,pr3.first,prob*pr1.second*pr2.second*pr3.second);
					};
				};
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

		// What ixns exist? for sampling
		_sampling_exists_w = false;
		_sampling_exists_h = false;
		_sampling_exists_j = false;
		_sampling_exists_k = false;

		// Does the lattice have the following structure?
		_latt_has_nn_structure = false;
		_latt_has_triplet_structure = false;
		_latt_has_quartic_structure = false;

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
		_sampling_exists_w = other._sampling_exists_w;
		_sampling_exists_h = other._sampling_exists_h;
		_sampling_exists_j = other._sampling_exists_j;
		_sampling_exists_k = other._sampling_exists_k;
		_latt_has_nn_structure = other._latt_has_nn_structure;
		_latt_has_triplet_structure = other._latt_has_triplet_structure;
		_latt_has_quartic_structure = other._latt_has_quartic_structure;
	};
	void Lattice::_reset() {
		_box_length = 0;
		_sp_map.clear();
		_sp_vec.clear();
		_latt.clear();
		_sampling_exists_w = false;
		_sampling_exists_h = false;
		_sampling_exists_j = false;
		_sampling_exists_k = false;
		_latt_has_nn_structure = false;
		_latt_has_triplet_structure = false;
		_latt_has_quartic_structure = false;
	};

	/********************
	Getters
	********************/

	int Lattice::dim() const {
		return _dim;
	};

	int Lattice::box_length() const {
		return _box_length;
	};

	/********************
	Add a species
	********************/

	void Lattice::add_species_possibility(Species *sp) {
		if (sp) { // not null
			_sp_map[sp->name()] = sp;
			_sp_vec.push_back(sp);

			// Add to all the sites
			for (latt_it lit=_latt.begin(); lit != _latt.end(); lit++) {
				lit->add_species_possibility(sp);
			};
		};
	};

	/********************
	Indicate that the hidden unit exists
	********************/

	void Lattice::set_sampling_flag_exists_w(bool flag) {
		_sampling_exists_w = flag;
	};
	void Lattice::set_sampling_flag_exists_h(bool flag) {
		_sampling_exists_h = flag;
	};
	void Lattice::set_sampling_flag_exists_j(bool flag) {
		_sampling_exists_j = flag;
	};
	void Lattice::set_sampling_flag_exists_k(bool flag) {
		if (_dim != 1) {
			std::cerr << "ERROR: triplets are only supported for d=1 currently" << std::endl;
			exit(EXIT_FAILURE);
		};
		_sampling_exists_k = flag;
	};

	/********************
	Initialize structure for NNs, triplets, etc.
	********************/

	void Lattice::init_nn_structure() {
		// Check: only do this once!
		if (!_latt_has_nn_structure) {

			// Set up neighbors
			std::vector<Site> nbrs;
			latt_it lit = _latt.begin();
			while (lit != _latt.end()) {
				// Neighbors
				nbrs.clear();
				if (lit->x-1 >= 1) {
					nbrs.push_back(Site(lit->x-1,lit->y,lit->z));
				};
				if (lit->x+1 <= _box_length) {
					nbrs.push_back(Site(lit->x+1,lit->y,lit->z));
				};
				if (_dim == 2 || _dim == 3) {
					if (lit->y-1 >= 1) {
						nbrs.push_back(Site(lit->x,lit->y-1,lit->z));
					};
					if (lit->y+1 <= _box_length) {
						nbrs.push_back(Site(lit->x,lit->y+1,lit->z));
					};
				};
				if (_dim == 3) {
					if (lit->z-1 >= 1) {
						nbrs.push_back(Site(lit->x,lit->y,lit->z-1));
					};
					if (lit->z+1 <= _box_length) {
						nbrs.push_back(Site(lit->x,lit->y,lit->z+1));
					};
				};

				// Go through neighbors
				for (auto nbr: nbrs) {
					// Add as nbr
					if (_dim == 1) {
						lit->nbrs.push_back(_look_up(nbr.x));
					} else if (_dim == 2) {
						lit->nbrs.push_back(_look_up(nbr.x,nbr.y));
					} else if (_dim == 3) {
						lit->nbrs.push_back(_look_up(nbr.x,nbr.y,nbr.z));
					};
				};

				// Next
				lit++;
			};

			// Now we have the structure
			_latt_has_nn_structure = true;
		};
	};

	void Lattice::init_triplet_structure() {
		// Check: only do this once!
		if (!_latt_has_triplet_structure) {

			latt_it lit1,lit2;
			for (latt_it lit = _latt.begin(); lit != _latt.end(); lit++) {
				if (lit->x != 1 && lit->x != 2) {
					// Both to the left
					lit1 = _look_up(lit->x-2);
					lit2 = _look_up(lit->x-1);
					lit->nbrs_triplets.push_back(LattIt2(lit1,lit2));
				};
				if (lit->x != 1 && lit->x != _box_length) {
					// One left, one right
					lit1 = _look_up(lit->x-1);
					lit2 = _look_up(lit->x+1);
					lit->nbrs_triplets.push_back(LattIt2(lit1,lit2));
				};
				if (lit->x != _box_length-1 && lit->x != _box_length) {
					// Both to the right
					lit1 = _look_up(lit->x+1);
					lit2 = _look_up(lit->x+2);
					lit->nbrs_triplets.push_back(LattIt2(lit1,lit2));
				};
			};

			// Now we have the structure
			_latt_has_triplet_structure = true;
		};
	};

	void Lattice::init_quartic_structure() {
		// Check: only do this once!
		if (!_latt_has_quartic_structure) {

			latt_it lit1,lit2,lit3;
			for (latt_it lit = _latt.begin(); lit != _latt.end(); lit++) {
				if (lit->x != 1 && lit->x != 2 && lit->x != 3) {
					// Three to the left
					lit1 = _look_up(lit->x-3);
					lit2 = _look_up(lit->x-2);
					lit3 = _look_up(lit->x-1);
					lit->nbrs_quartics.push_back(LattIt3(lit1,lit2,lit3));
				};
				if (lit->x != 1 && lit->x != 2 && lit->x != _box_length) {
					// Two to the left, one to the right
					lit1 = _look_up(lit->x-2);
					lit2 = _look_up(lit->x-1);
					lit3 = _look_up(lit->x+1);
					lit->nbrs_quartics.push_back(LattIt3(lit1,lit2,lit3));
				};
				if (lit->x != 1 && lit->x != _box_length-1 && lit->x != _box_length) {
					// One left, two to the right
					lit1 = _look_up(lit->x-1);
					lit2 = _look_up(lit->x+1);
					lit3 = _look_up(lit->x+2);
					lit->nbrs_quartics.push_back(LattIt3(lit1,lit2,lit3));
				};
				if (lit->x != _box_length-2 && lit->x != _box_length-1 && lit->x != _box_length) {
					// Three to the right
					lit1 = _look_up(lit->x+1);
					lit2 = _look_up(lit->x+2);
					lit3 = _look_up(lit->x+3);
					lit->nbrs_quartics.push_back(LattIt3(lit1,lit2,lit3));
				};
			};

			// Now we have the structure
			_latt_has_quartic_structure = true;
		};
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
			it->set_site_empty();
		};

		// Reset the counts and nns
		for (auto sp: _sp_vec) {
			sp->reset_counts();
		};
	};
	int Lattice::size() { return _latt.size(); };

	/********************
	Binarize
	********************/

	void Lattice::binarize() {
		for (latt_it lit=_latt.begin(); lit != _latt.end(); lit++) {
			lit->binarize();
		};
	};

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

	void Lattice::read_from_file(std::string fname, bool binary)
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
		std::string prob="";
		double prob_val;
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
			    if (!binary) {
			    	// Read the prob
			    	iss >> prob;
			    };
		    	// Add to lattice
		    	if (_dim == 1) {
		    		lit = _look_up(atoi(x.c_str()));
			    } else if (_dim == 2) {
		    		lit = _look_up(atoi(x.c_str()),atoi(y.c_str()));
			    } else if (_dim == 3) {
			    	lit = _look_up(atoi(x.c_str()),atoi(y.c_str()),atoi(z.c_str()));
			    };
			    if (binary) {
		    		lit->set_prob(_sp_map[sp],1.0);
		    	} else {
		    		prob_val = atof(prob.c_str());
		    		lit->set_prob(_sp_map[sp],prob_val);
		    		lit->set_prob(nullptr,1.0-prob_val);
		    	};
	    		// Reset
		    	sp=""; x=""; y=""; z=""; prob="";
			};
		};
		f.close();

		// std::cout << "Read: " << fname << std::endl;
	};

	/********************
	Populate randomly according to some counts
	********************/

	void Lattice::populate_randomly(std::map<Species*, int> counts) {
		// Clear the current lattice
		clear();

		bool did_place;
		int ctr_tries;
		latt_it lit;

		// Go through the species
		for (auto pr: counts) {
			// Go through counts
			for (int i=0; i<pr.second; i++) {
				// Try to place
				did_place = false;
				ctr_tries = 0;
				while (did_place == false && ctr_tries < 1000) { // Try 1000 different places
			    	if (_dim == 1) {
			    		lit = _look_up(randI(1,_box_length));
				    } else if (_dim == 2) {
			    		lit = _look_up(randI(1,_box_length),randI(1,_box_length));
				    } else if (_dim == 3) {
			    		lit = _look_up(randI(1,_box_length),randI(1,_box_length),randI(1,_box_length));
				    };
				    // Check if empty
				    if (lit->empty()) {
				    	// Yes, it's empty - place!
				    	did_place = true;
				    	lit->set_prob(pr.first,1.0);
				    } else {
				    	// Try again!
				    	ctr_tries++;
				    };
				};
				// Check we didn't run out
				if (ctr_tries >= 1000) {
					std::cerr << "WARNING! Couldn't place all the mols for species: " << pr.first->name() << " wanted to place: " << pr.second << " mols but only got to: " << i << std::endl;
					// Don't try to place any more
					break;
				};
			};
		};
	};

	/********************
	Sample
	********************/

	void Lattice::sample(bool binary) {

		// Declarations
		latt_it it;
		double energy;
		std::map<Species*, std::vector< std::pair< HiddenUnit*, std::vector<IxnParam*> > > >::iterator it_hups;
		int i_chosen;
		std::vector<double> props,probs;
		latt_it lit;
		double tot;

		// Go through all lattice sites
		for (it = _latt.begin(); it != _latt.end(); it++) {

			// Clear props, probs
			props.clear();
			probs.clear();
			props.push_back(0.0);

			// Empty = 1
			props.push_back(1.0);
			probs.push_back(1.0);

			// Go through all possible species this could be, calculate propensities
			for (auto sp_new: _sp_vec) {
				
				// Get energy
				energy = it->get_activation(sp_new);

				// Append prop
				props.push_back(props.back()+exp(energy));
				probs.push_back(exp(energy));
			};

			/*
			for (auto pr: props) {
				std::cout << pr << " ";
			};
			std::cout << std::endl;
			*/

			// Commit probs
			if (binary) {

				// Sample RV
				// i_chosen = _sample_prob_vec(probs);
				i_chosen = sample_prop_vec(props);

				if (i_chosen==0) {
					// Flip down (new spin = 0)
					it->set_site_empty();
				} else {
					// Make the appropriate species at this site (guaranteed empty)
					it->set_site_binary(_sp_vec[i_chosen-1]);
				};

			} else {

				// Normalize probs
				tot=0.0;
				for (auto pr: probs) {
					tot += pr;
				};

				// Write into species
				it->set_prob(nullptr,probs[0]/tot);
				for (int i=0; i<_sp_vec.size(); i++) {
					it->set_prob(_sp_vec[i],probs[i+1]/tot);
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
};