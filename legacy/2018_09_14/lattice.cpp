#include "../../include/dynamicboltz_bits/Latt.hpp"

// Other headers
#include "../../include/dynamicboltz_bits/general.hpp"
#include "site.hpp"
#include "species.hpp"
#include "ixn_param_traj.hpp"

#include <iostream>
#include <fstream>
#include <numeric>
#include "math.h"
#include <ctime>
#include <sstream>
#include <random>

/************************************
* Namespace for dblz
************************************/

namespace dblz {

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

		// Does the Latt have the following structure?
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
	Validate graph
	********************/

	void Lattice::validate_graph() const {
		// Go through all sites
		for (auto const &s: _latt) {
			s.validate();
		};
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
			for (Latt_it lit=_latt.begin(); lit != _latt.end(); lit++) {
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
			Latt_it lit = _latt.begin();
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
			for (Latt_it lit = _latt.begin(); lit != _latt.end(); lit++) {
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
			for (Latt_it lit = _latt.begin(); lit != _latt.end(); lit++) {
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
		// Clear the Latt
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
		for (Latt_it lit=_latt.begin(); lit != _latt.end(); lit++) {
			lit->binarize();
		};
	};

	/********************
	Write Latt to a file
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
	Read Latt from a file
	********************/

	void Lattice::read_from_file(std::string fname, bool binary)
	{
		// Clear the current Latt
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
		    	// Add to Latt
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
		// At most half the Latt is filled
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
		// Clear the current Latt
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
		Latt_it it = _latt.begin();
		std::advance(it,n);
		return &*it;
	};
	Site* Lattice::_look_up(int x, int y) {
		// Figure out index in list
		int n = (x-1)*_box_length + y-1;

		// Grab
		Latt_it it = _latt.begin();
		std::advance(it,n);
		return &*it;
	};
	Site* Lattice::_look_up(int x, int y, int z) {
		// Figure out index in list
		int n = (x-1)*_box_length*_box_length + (y-1)*_box_length + z-1;

		// Grab
		Latt_it it = _latt.begin();
		std::advance(it,n);
		return &*it;
	};
};