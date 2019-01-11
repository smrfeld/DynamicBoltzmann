#include "lattice.hpp"
#include <iostream>
#include <fstream>
#include <numeric>
#include "general.hpp"
#include "math.h"
#include <ctime>
#include "species.hpp"

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Struct to hold a lattice Site
	****************************************/

	// Constructor
	Site::Site() : Site(0,0,0,nullptr) {};
	Site::Site(int xIn, int yIn, int zIn) : Site(xIn, yIn, zIn, nullptr) {};	
	Site::Site(int xIn, int yIn, int zIn, Species *spIn) {
		this->x = xIn;
		this->y = yIn;
		this->z = zIn;
		this->sp = spIn;
	};	

	// Comparator
	bool operator <(const Site& a, const Site& b) {
    	return std::tie(a.x, a.y, a.z) < std::tie(b.x, b.y, b.z);
	};
	bool operator==(const Site& a, const Site& b) {
		return std::tie(a.x, a.y, a.z) == std::tie(b.x, b.y, b.z);
	}; 
	std::ostream& operator<<(std::ostream& os, const Site& s)
	{
		if (s.sp) {
		    return os << s.x << " " << s.y << " " << s.z << " " << s.sp->name;
		} else {
		    return os << s.x << " " << s.y << " " << s.z << " empty";
		};
	};

	/****************************************
	Lattice
	****************************************/

	/********************
	Constructor/Destructor
	********************/

	// Constructor
	Lattice::Lattice(int box_length)
	{
		this->_box_length = box_length;

		// Make a fully linked list of sites
		for (int x=1; x<=box_length; x++) {
			for (int y=1; y<=box_length; y++) {
				for (int z=1; z<=box_length; z++) {
					_latt.push_back(Site(x,y,z));					
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
			if (lit->y-1 >= 1) {
				nbrs.push_back(Site(lit->x,lit->y-1,lit->z));
			};
			if (lit->y+1 <= box_length) {
				nbrs.push_back(Site(lit->x,lit->y+1,lit->z));
			};
			if (lit->z-1 >= 1) {
				nbrs.push_back(Site(lit->x,lit->y,lit->z-1));
			};
			if (lit->z+1 <= box_length) {
				nbrs.push_back(Site(lit->x,lit->y,lit->z+1));
			};

			// Go through neighbors
			for (auto nbr: nbrs) {
				// Add as nbr
				lit->nbrs.push_back(_look_up(nbr.x,nbr.y,nbr.z));
			};

			// Next
			lit++;
		};
	};
	Lattice::Lattice()
	{
		this->_box_length = 0;
	};

	// Destructor
	Lattice::~Lattice() {};

	/********************
	Add a species
	********************/

	void Lattice::add_species(Species *sp) {
		if (sp) {
			_sp_map[sp->name] = sp;
		};
	};

	/********************
	Clear, size
	********************/

	void Lattice::clear() { 
		// Reset the counts and nns
		for (auto sp_pair: _sp_map) {
			sp_pair.second->count = 0;
			for (auto sp_pair2: _sp_map) {
				sp_pair.second->nn_count[sp_pair2.second] = 0;
			};
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
		std::cout << "Making: " << s->x << " " << s->y << " " << s->z;
		int ctr=0;
		for (auto nbr_it: s->nbrs) { // go through nbrs
			if (nbr_it->sp) { // is nbr site occ
				ctr++;
			};
		};
		std::cout << " no nbrs: " << ctr << " count now: " << sp->nn_count[_sp_map["A"]];
		*/

		// Make
		s->sp = sp;

		// Update count and nns
		sp->count++;
		for (auto nbr_it: s->nbrs) { // go through nbrs
			if (nbr_it->sp) { // is nbr site occ
				// Update bidirectionally, unless it's the same species
				sp->nn_count[nbr_it->sp]++;
				if (nbr_it->sp != sp) {
					nbr_it->sp->nn_count[sp]++;
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
		s->sp->count--;
		for (auto nbr_it: s->nbrs) { // go through nbrs
			if (nbr_it->sp) { // is nbr site occ
				// Update bidirectionally, unless it's the same species
				s->sp->nn_count[nbr_it->sp]--;
				if (nbr_it->sp != s->sp) {
					nbr_it->sp->nn_count[s->sp]--;
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
				f << l.x << " " << l.y << " " << l.z << " " << l.sp->name << "\n";
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
		char frag[100]; // fragments of the line
		int i_frag=0;
		std::string x="",y="",z="";
		std::string sp="";
		if (f.is_open()) { // make sure we found it
			while (!f.eof()) {
			    f >> frag;
			    if (i_frag==0) {
			    	x += frag; i_frag++;
			    } else if (i_frag==1) {
			    	y += frag; i_frag++;
			    } else if (i_frag==2) {
			    	z += frag; i_frag++;
			    } else if (i_frag==3) {
			    	sp += frag;
			    	// Add to map
			    	make_mol(_look_up(atoi(x.c_str()),atoi(y.c_str()),atoi(z.c_str())),_sp_map[sp]);
			    	i_frag = 0;
			    	sp=""; x=""; y=""; z="";
			    };
			};
		};
		f.close();

		// std::cout << "Read: " << fname << " no nn: " << _sp_map["A"]->nn_count[_sp_map["A"]] << std::endl;
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

	latt_it Lattice::_look_up(int x, int y, int z) {
		// Figure out index in list
		int n = (x-1)*_box_length*_box_length + (y-1)*_box_length + z-1;

		// Grab
		latt_it it = _latt.begin();
		std::advance(it,n);
		return it;
	};

};