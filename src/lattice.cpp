#include "ixn_param_traj.hpp" // also includes lattice header
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
	ConnectionVH::ConnectionVH(Site *site, HiddenUnit *hidden_unit, std::vector<IxnParamTraj*> ips) {
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
	void ConnectionVH::add_ixn_param(IxnParamTraj* ip) {
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
	Class to hold a lattice Site
	****************************************/

	/********************
	Constructor
	********************/

	Site::Site(int x) : Site(x,0,0) { _dim=1; };
	Site::Site(int x, int y) : Site(x,y,0) { _dim=2; };
	Site::Site(int x, int y, int z) {
		_dim=3;
		_x = x;
		_y = y;
		_z = z;
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
		_dim = 0;
		_x = 0;
		_y = 0;
		_z = 0;
		_nbrs.clear();
		_nbrs_triplets.clear();
		_nbrs_quartics.clear();
		_hidden_conns.clear();
		_prob_empty = 0.0;
		_probs.clear();
		_sp_possible.clear();
	};
	void Site::_copy(const Site& other) {
		_dim = other._dim;
		_x = other._x;
		_y = other._y;
		_z = other._z;
		_nbrs = other._nbrs;
		_nbrs_triplets = other._nbrs_triplets;
		_nbrs_quartics = other._nbrs_quartics;
		_hidden_conns = other._hidden_conns;
		_prob_empty = other._prob_empty;
		_probs = other._probs;
		_sp_possible = other._sp_possible;
	};

	/********************
	Check location
	********************/

	int Site::x() const {
		return _x;
	};
	int Site::y() const {
		return _y;
	};
	int Site::z() const {
		return _z;
	};
	bool Site::less_than(const Site &other) const {
		if (_dim == 1) {
			return _x < other._x;
		} else if (_dim == 2) {
			return std::tie(_x, _y) < std::tie(other._x, other._y);
 		} else if (_dim == 3) {
			return std::tie(_x, _y, _z) < std::tie(other._x, other._y, other._z);
	    } else {
	    	return false;
	    };	
	};

	/********************
	Add neighbors
	********************/

	void Site::add_nbr(Site *s) {
		_nbrs.push_back(s);
	};
	void Site::add_nbr_triplet(Site *s1, Site *s2) {
		_nbrs_triplets.push_back(Site2(s1,s2));
	};
	void Site::add_nbr_quartic(Site *s1, Site *s2, Site *s3) {
		_nbrs_quartics.push_back(Site3(s1,s2,s3));
	};

	/********************
	Comparator
	********************/

	bool operator <(const Site& a, const Site& b) {
		return a.less_than(b);
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
		_sp_possible.push_back(sp);
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
	Sample
	********************/

	void Site::sample(bool binary) {
		double energy;

		std::vector<double> props,probs;
		props.push_back(0.0);

		// Empty = 1
		props.push_back(1.0);
		probs.push_back(1.0);

		// Go through all possible species this could be, calculate propensities
		for (auto sp: _sp_possible) {
			
			// Bias
			energy = sp->h();

			// NNs - go through neihbors
			for (auto lit: _nbrs) {
				// Get all probs
				const std::map<Species*, double> prs = lit->get_probs();
				// Go through all probs
				for (auto pr: prs) {
					// J * prob
					energy += sp->j(pr.first) * pr.second;
				};
			};

			// Triplets - go through all pairs to consider
			for (auto trip: _nbrs_triplets) {
				// Get all probs
				const std::map<Species*, double> prs1 = trip.s1->get_probs();
				const std::map<Species*, double> prs2 = trip.s2->get_probs();
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

			// Append prop
			energy = exp(energy);
			props.push_back(props.back()+energy);
			probs.push_back(energy);
		};

		// Commit probs
		if (binary) {

			// Sample RV
			int i_chosen = sample_prop_vec(props);

			if (i_chosen==0) {
				// Flip down (new spin = 0)
				set_site_empty();
			} else {
				// Make the appropriate species at this site (guaranteed empty)
				set_site_binary(_sp_possible[i_chosen-1]);
			};

		} else {

			// Normalize probs
			double tot=0.0;
			for (auto pr: probs) {
				tot += pr;
			};

			// Write into species
			set_prob(nullptr,probs[0]/tot);
			for (int i=0; i<_sp_possible.size(); i++) {
				set_prob(_sp_possible[i],probs[i+1]/tot);
			};
		};

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
		for (auto nbr_it: _nbrs) {
			// Get all the probs
			const std::map<Species*, double> prs = nbr_it->get_probs();
			for (auto pr: prs) {
				// Increment
				sp->nn_count_increment(pr.first,prob*pr.second);
			};
		};

		// Triplets, if needed
		for (auto trip: _nbrs_triplets) {
			// Get all the probs
			const std::map<Species*, double> prs1 = trip.s1->get_probs();
			const std::map<Species*, double> prs2 = trip.s2->get_probs();
			for (auto pr1: prs1) {
				for (auto pr2: prs2) {
					// Increment
					sp->triplet_count_increment(pr1.first,pr2.first,prob*pr1.second*pr2.second);
				};
			};
		};

		// Quartics, if needed
		for (auto quart: _nbrs_quartics) {
			// Get all the probs
			const std::map<Species*, double> prs1 = quart.s1->get_probs();
			const std::map<Species*, double> prs2 = quart.s2->get_probs();
			const std::map<Species*, double> prs3 = quart.s3->get_probs();
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
		_latt_has_nn_structure = other._latt_has_nn_structure;
		_latt_has_triplet_structure = other._latt_has_triplet_structure;
		_latt_has_quartic_structure = other._latt_has_quartic_structure;
	};
	void Lattice::_reset() {
		_box_length = 0;
		_sp_map.clear();
		_sp_vec.clear();
		_latt.clear();
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
				if (lit->x()-1 >= 1) {
					nbrs.push_back(Site(lit->x()-1,lit->y(),lit->z()));
				};
				if (lit->x()+1 <= _box_length) {
					nbrs.push_back(Site(lit->x()+1,lit->y(),lit->z()));
				};
				if (_dim == 2 || _dim == 3) {
					if (lit->y()-1 >= 1) {
						nbrs.push_back(Site(lit->x(),lit->y()-1,lit->z()));
					};
					if (lit->y()+1 <= _box_length) {
						nbrs.push_back(Site(lit->x(),lit->y()+1,lit->z()));
					};
				};
				if (_dim == 3) {
					if (lit->z()-1 >= 1) {
						nbrs.push_back(Site(lit->x(),lit->y(),lit->z()-1));
					};
					if (lit->z()+1 <= _box_length) {
						nbrs.push_back(Site(lit->x(),lit->y(),lit->z()+1));
					};
				};

				// Go through neighbors
				for (auto nbr: nbrs) {
					// Add as nbr
					if (_dim == 1) {
						lit->add_nbr(_look_up(nbr.x()));
					} else if (_dim == 2) {
						lit->add_nbr(_look_up(nbr.x(),nbr.y()));
					} else if (_dim == 3) {
						lit->add_nbr(_look_up(nbr.x(),nbr.y(),nbr.z()));
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

			Site *s1,*s2;
			for (latt_it lit = _latt.begin(); lit != _latt.end(); lit++) {
				if (lit->x() != 1 && lit->x() != 2) {
					// Both to the left
					s1 = _look_up(lit->x()-2);
					s2 = _look_up(lit->x()-1);
					lit->add_nbr_triplet(s1,s2);
				};
				if (lit->x() != 1 && lit->x() != _box_length) {
					// One left, one right
					s1 = _look_up(lit->x()-1);
					s2 = _look_up(lit->x()+1);
					lit->add_nbr_triplet(s1,s2);
				};
				if (lit->x() != _box_length-1 && lit->x() != _box_length) {
					// Both to the right
					s1 = _look_up(lit->x()+1);
					s2 = _look_up(lit->x()+2);
					lit->add_nbr_triplet(s1,s2);
				};
			};

			// Now we have the structure
			_latt_has_triplet_structure = true;
		};
	};

	void Lattice::init_quartic_structure() {
		// Check: only do this once!
		if (!_latt_has_quartic_structure) {

			Site *s1,*s2,*s3;
			for (latt_it lit = _latt.begin(); lit != _latt.end(); lit++) {
				if (lit->x() != 1 && lit->x() != 2 && lit->x() != 3) {
					// Three to the left
					s1 = _look_up(lit->x()-3);
					s2 = _look_up(lit->x()-2);
					s3 = _look_up(lit->x()-1);
					lit->add_nbr_quartic(s1,s2,s3);
				};
				if (lit->x() != 1 && lit->x() != 2 && lit->x() != _box_length) {
					// Two to the left, one to the right
					s1 = _look_up(lit->x()-2);
					s2 = _look_up(lit->x()-1);
					s3 = _look_up(lit->x()+1);
					lit->add_nbr_quartic(s1,s2,s3);
				};
				if (lit->x() != 1 && lit->x() != _box_length-1 && lit->x() != _box_length) {
					// One left, two to the right
					s1 = _look_up(lit->x()-1);
					s2 = _look_up(lit->x()+1);
					s3 = _look_up(lit->x()+2);
					lit->add_nbr_quartic(s1,s2,s3);
				};
				if (lit->x() != _box_length-2 && lit->x() != _box_length-1 && lit->x() != _box_length) {
					// Three to the right
					s1 = _look_up(lit->x()+1);
					s2 = _look_up(lit->x()+2);
					s3 = _look_up(lit->x()+3);
					lit->add_nbr_quartic(s1,s2,s3);
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
		return _look_up(x);
	};
	Site* Lattice::get_site(int x, int y) {
		if (_dim != 2) {
			std::cerr << "ERROR: dim wrong in get_site" << std::endl;
			exit(EXIT_FAILURE);
		};
		return _look_up(x,y);
	};
	Site* Lattice::get_site(int x, int y, int z) {
		if (_dim != 3) {
			std::cerr << "ERROR: dim wrong in get_site" << std::endl;
			exit(EXIT_FAILURE);
		};
		return _look_up(x,y,z);
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
			if (_dim == 1) {
				f << l.x() << "\n";
			} else if (_dim == 2) {
				f << l.x() << " " << l.y() << "\n";
			} else if (_dim == 3) {
				f << l.x() << " " << l.y() << " " << l.z() << "\n";
			};
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
		Site* s;
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
		    		s = _look_up(atoi(x.c_str()));
			    } else if (_dim == 2) {
		    		s = _look_up(atoi(x.c_str()),atoi(y.c_str()));
			    } else if (_dim == 3) {
			    	s = _look_up(atoi(x.c_str()),atoi(y.c_str()),atoi(z.c_str()));
			    };
			    if (binary) {
		    		s->set_prob(_sp_map[sp],1.0);
		    	} else {
		    		prob_val = atof(prob.c_str());
		    		s->set_prob(_sp_map[sp],prob_val);
		    		s->set_prob(nullptr,1.0-prob_val);
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

	void Lattice::populate_randomly() {
		// Random number of initial particles (min is 1, max is box vol)
		int n = randI(1, pow(_box_length,_dim));

		// Random initial counts
		// Don't populate too much, else this is hard to find empty sites to place mols!
		// At most half the lattice is filled
		int n_possible = pow(_box_length,_dim) / 2;
		std::map<Species*,int> counts;
		for (auto sp: _sp_vec) {
			counts[sp] = randI(0,n_possible);
			n_possible -= counts[sp];
			if (n_possible < 0) { n_possible = 0; };
		};

		// Populate at random positions
		populate_randomly(counts);
	};

	void Lattice::populate_randomly(std::map<Species*, int> counts) {
		// Clear the current lattice
		clear();

		bool did_place;
		int ctr_tries;
		Site *s;

		// Go through the species
		for (auto pr: counts) {
			// Go through counts
			for (int i=0; i<pr.second; i++) {
				// Try to place
				did_place = false;
				ctr_tries = 0;
				while (did_place == false && ctr_tries < 1000) { // Try 1000 different places
			    	if (_dim == 1) {
			    		s = _look_up(randI(1,_box_length));
				    } else if (_dim == 2) {
			    		s = _look_up(randI(1,_box_length),randI(1,_box_length));
				    } else if (_dim == 3) {
			    		s = _look_up(randI(1,_box_length),randI(1,_box_length),randI(1,_box_length));
				    };
				    // Check if empty
				    if (s->empty()) {
				    	// Yes, it's empty - place!
				    	did_place = true;
				    	s->set_prob(pr.first,1.0);
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

		for (auto it = _latt.begin(); it != _latt.end(); it++) {

			it->sample(binary);

		};
	};

	/****************************************
	Lattice - PRIVATE
	****************************************/

	/********************
	Lookup a site iterator from x,y,z
	********************/

	Site* Lattice::_look_up(int x) {
		// Figure out index in list
		int n = x-1;

		// Grab
		latt_it it = _latt.begin();
		std::advance(it,n);
		return &*it;
	};
	Site* Lattice::_look_up(int x, int y) {
		// Figure out index in list
		int n = (x-1)*_box_length + y-1;

		// Grab
		latt_it it = _latt.begin();
		std::advance(it,n);
		return &*it;
	};
	Site* Lattice::_look_up(int x, int y, int z) {
		// Figure out index in list
		int n = (x-1)*_box_length*_box_length + (y-1)*_box_length + z-1;

		// Grab
		latt_it it = _latt.begin();
		std::advance(it,n);
		return &*it;
	};
};