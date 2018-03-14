#include "lattice.hpp"
#include <iostream>
#include <fstream>
#include <numeric>
#include "general.hpp"
#include "math.h"
#include <ctime>
#include "species.hpp"
#include <sstream>

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
	};	
	Site::Site(const Site& other) {
		dim = other.dim;
		x = other.x;
		y = other.y;
		z = other.z;
		sp = other.sp;
		nbrs = other.nbrs;
	};
	Site::Site(Site&& other) {
		dim = other.dim;
		x = other.x;
		y = other.y;
		z = other.z;
		sp = other.sp;
		nbrs = other.nbrs;
		// Clear other
		other.dim = 0;
		other.x = 0;
		other.y = 0;
		other.z = 0;
		other.sp = nullptr;
		other.nbrs.clear();
	};
	Site& Site::operator=(const Site& other) {
		if (this != &other) {
			dim = other.dim;
			x = other.x;
			y = other.y;
			z = other.z;
			sp = other.sp;
			nbrs = other.nbrs;
		};
		return *this;
	};
	Site& Site::operator=(Site&& other) {
		if (this != &other) {
			dim = other.dim;
			x = other.x;
			y = other.y;
			z = other.z;
			sp = other.sp;
			nbrs = other.nbrs;
			// Clear other
			other.dim = 0;
			other.x = 0;
			other.y = 0;
			other.z = 0;
			other.sp = nullptr;
			other.nbrs.clear();
		};
		return *this;
	};
	Site::~Site() {
		// Nothing...
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
	};
	void Lattice::_copy(Lattice &&other) {
		_latt = other._latt;
		_box_length = other._box_length;
		_sp_map = other._sp_map;
		// Clear other
		other._box_length = 0;
		other._sp_map.clear();
		other._latt.clear();
	};

	/********************
	Add a species
	********************/

	void Lattice::add_species(Species *sp) {
		if (sp) {
			_sp_map[sp->name()] = sp;
		};
	};

	/********************
	Clear, size
	********************/

	void Lattice::clear() { 
		// Reset the counts and nns
		for (auto sp_pair: _sp_map) {
			sp_pair.second->reset_counts();
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
		sp->count_plus();
		for (auto nbr_it: s->nbrs) { // go through nbrs
			if (nbr_it->sp) { // is nbr site occ
				// Update bidirectionally, unless it's the same species
				sp->nn_count_plus(nbr_it->sp);
				if (nbr_it->sp != sp) {
					nbr_it->sp->nn_count_plus(sp);
				};
			};
		};

		// std::cout << " now: count: " << sp->nn_count[_sp_map["A"]] << std::endl;

		return true;
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
		s->sp->count_minus();
		for (auto nbr_it: s->nbrs) { // go through nbrs
			if (nbr_it->sp) { // is nbr site occ
				// Update bidirectionally, unless it's the same species
				s->sp->nn_count_minus(nbr_it->sp);
				if (nbr_it->sp != s->sp) {
					nbr_it->sp->nn_count_minus(s->sp);
				};
			};
		};

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
	Anneal
	********************/

	void Lattice::anneal(int n_steps) {

		// Annealing temperature
		// double temp0 = 3.0;
		// double temp;
		// Formula:
		// temp0 * exp( - log(temp0) * i / (n_steps-1) )
		// Such that final temp is 1

		// Timing
		
		// std::clock_t    start;

		// Declarations

		double hOld,jOld,hNew,jNew,energy_diff;
		latt_it it_flip;
		std::map<std::string,Species*>::iterator it_sp;
		Species *sp_new = nullptr;

		// Go through the steps
		for (int i=0; i<n_steps; i++) {

			// Annealing temp
			// temp = temp0 * exp( - log(temp0) * i / (n_steps-1) );

			// start = std::clock();

			// Pick a site to flip randomly
			it_flip = _latt.begin();
			std::advance(it_flip,randI(0,_latt.size()-1));

			//std::cout << "NBRS Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
			//start = std::clock();

			// Check if this site is occupied
			if (it_flip->sp) {
				// Occupied - flip down
				hOld = -it_flip->sp->h();
				jOld = 0.0;
				for (auto it_nbr: it_flip->nbrs) {
					if (it_nbr->sp) {
						jOld -= it_flip->sp->j(it_nbr->sp);
					};
				};
				// New couplings
				hNew = 0.0;
				jNew = 0.0;
			} else {
				// Unoccupied - flip up
				hOld = 0.0;
				jOld = 0.0;
				// Random species
				it_sp = _sp_map.begin();
				std::advance(it_sp, randI(0,_sp_map.size()-1));
				sp_new = it_sp->second;
				// New couplings
				hNew = -sp_new->h();
				jNew = 0.0;
				for (auto it_nbr: it_flip->nbrs) {
					if (it_nbr->sp) {
						jNew -= sp_new->j(it_nbr->sp);
					};
				};
			};

			//std::cout << "CALC Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
			//start = std::clock();

			// Energy difference
			energy_diff = hNew + jNew - hOld - jOld;
			if (energy_diff < 0.0 || exp(-energy_diff) > randD(0.0,1.0)) {
				// Accept the flip!
				if (it_flip->sp) {
					// Occupied - flip down
					erase_mol(it_flip);
				} else {
					// Unoccupied - flip up
					make_mol(it_flip,sp_new);
				};
			};

			//std::cout << "FLIP Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
			//start = std::clock();
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