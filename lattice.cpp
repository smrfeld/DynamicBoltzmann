#include "lattice.hpp"
#include <iostream>
#include <fstream>
#include <numeric>
#include "general.hpp"
#include "math.h"

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Ordered string pair
	****************************************/

	StrPair::StrPair(std::string r1, std::string r2)
	{
		if (r1 < r2) {
			s1 = r1;
			s2 = r2;
		} else {
			s1 = r2;
			s2 = r1;
		};
	};
	bool operator <(const StrPair& a, const StrPair& b) {
    	return std::tie(a.s1, a.s2) < std::tie(b.s1, b.s2);
	};
	bool operator==(const StrPair& a, const StrPair& b) {
		return std::tie(a.s1, a.s2) == std::tie(b.s1, b.s2);
	};

	/****************************************
	Structure to hold a lattice site iterator
	****************************************/

	SiteIt::SiteIt() {};
	SiteIt::SiteIt(lattice_map::iterator itIn, lattice_map_1::iterator it_1In, lattice_map_2::iterator it_2In) 
	{
		this->it = itIn;
		this->it_1 = it_1In;
		this->it_2 = it_2In;
	};
	std::ostream& operator<<(std::ostream& os, const SiteIt& sit)
	{
	    return os << sit.it->first << " " << sit.it_1->first << " " << sit.it_2->first;
	};

	/****************************************
	Struct to hold a lattice Site
	****************************************/

	// Constructor
	Site::Site() {};
	Site::Site(int xIn, int yIn, int zIn) {
		this->x = xIn;
		this->y = yIn;
		this->z = zIn;
	};
	Site::Site(SiteIt sit) {
		this->x = sit.it->first;
		this->y = sit.it_1->first;
		this->z = sit.it_2->first;
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
	    return os << s.x << " " << s.y << " " << s.z;
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
	};
	Lattice::Lattice()
	{
		this->_box_length = 0;
	};

	// Destructor
	Lattice::~Lattice() {};

	/********************
	Clear, size
	********************/

	void Lattice::clear() { this->_map.clear(); };
	int Lattice::size() { return this->_map.size(); };

	/********************
	Make a mol
	********************/

	std::pair<bool,SiteIt> Lattice::make_mol(Site s, std::string sp) 
	{
		// Check if the site is empty
		std::pair<bool,SiteIt> spair = get_mol_it(s);
		if (spair.first) {
			// Not empty
			return std::make_pair(false,SiteIt());
		};

		// Update counts
		// To-do?
		// sp->count++;

		// Make
		lattice_map::iterator it;
		lattice_map_1::iterator it_1;
		lattice_map_2::iterator it_2;
		it = this->_map.find(s.x);
		if (it == this->_map.end()) {
			auto ret = this->_map.insert(std::make_pair(s.x,lattice_map_1()));
			it = ret.first;
		};
		it_1 = it->second.find(s.y);
		if (it_1 == it->second.end()) {
			auto ret_1 = it->second.insert(std::make_pair(s.y,lattice_map_2()));
			it_1 = ret_1.first;
		};
		auto ret_2 = it_1->second.insert(std::make_pair(s.z,sp));
		it_2 = ret_2.first;

		return std::make_pair(true,SiteIt(it,it_1,it_2));
	};

	std::pair<bool,SiteIt> Lattice::make_mol_random(std::string sp) 
	{
		// Get a random free site
		std::pair<bool,Site> spair = get_free_site();

		if (!(spair.first)) {
			// No free sites at all
			return std::make_pair(false,SiteIt());
		};

		// Make
		std::pair<bool,SiteIt> ret = make_mol(spair.second,sp);
		return ret;
	};

	/********************
	Erase a mol
	********************/

	bool Lattice::erase_mol(Site s) 
	{
		// Get an iterator to the site
		std::pair<bool,SiteIt> spair = get_mol_it(s);
		if (spair.first) {
			return erase_mol_it(spair.second);
		} else {
			return false;
		};
	};

	bool Lattice::erase_mol_it(SiteIt sit) 
	{
		// Update counts on species
		// To-do?
		// sit.it_2->second.sp->count--;

		// Erase inner_2 from inner_1
		sit.it_1->second.erase(sit.it_2);
		if (sit.it_1->second.size() == 0)
		{
			// Erase inner_1 from outer
			sit.it->second.erase(sit.it_1);
			if (sit.it->second.size() == 0)
			{
				// Erase outer from map
				this->_map.erase(sit.it);
				return true;
			};
		};
		return true; // Nonsense; it can't fail?!
	};

	std::pair<bool,Site> Lattice::erase_mol_random(std::string sp) 
	{
		// Get an iterator to a random site
		std::pair<bool,SiteIt> spair = get_mol_random_it(sp);
		if (spair.first) {
			Site s = Site(spair.second.it->first,spair.second.it_1->first,spair.second.it_2->first);
			bool succ = erase_mol_it(spair.second);
			if (succ) {
				return std::make_pair(true,s);
			};
		};
		return std::make_pair(false,Site());
	};

	/********************
	Get a mol
	********************/

	std::pair<bool,SiteIt> Lattice::get_mol_it(Site s) 
	{
		auto it = this->_map.find(s.x);
		if (it != this->_map.end()) {
			auto it_1 = it->second.find(s.y);
			if (it_1 != it->second.end()) {
				auto it_2 = it_1->second.find(s.z);
				if (it_2 != it_1->second.end()) {
					return std::make_pair(true,SiteIt(it,it_1,it_2));
				};
			};
		};
		return std::make_pair(false,SiteIt());	
	};

	std::pair<bool,SiteIt> Lattice::get_mol_it(Site s, std::string sp) {
		std::pair<bool,SiteIt> ret = get_mol_it(s);
		if (ret.first) {
			if (ret.second.it_2->second == sp) {
				return ret;
			};
		};
		return std::make_pair(false,SiteIt());
	};

	std::pair<bool,SiteIt> Lattice::get_mol_random_it() {
		// Get random indexes to search
		std::map<int,std::vector<int>> idxs = _get_random_idxs();

		// Try all sites
		std::pair<bool,SiteIt> spair;
		for (auto i1=0; i1 < this->_box_length; i1++) {
			for (auto i2=0; i2 < this->_box_length; i2++) {
				for (auto i3=0; i3 < this->_box_length; i3++) {
					spair = get_mol_it(Site(idxs[0][i1],idxs[1][i2],idxs[2][i3]));
					if (spair.first) 
					{
						return spair;
					};
				};
			};
		};
		return std::make_pair(false,SiteIt());
	};

	std::pair<bool,SiteIt> Lattice::get_mol_random_it(std::string sp) 
	{
		// Get random indexes to search
		std::map<int,std::vector<int>> idxs = _get_random_idxs();

		// Try all sites
		std::pair<bool,SiteIt> spair;
		for (auto i1=0; i1 < this->_box_length; i1++) {
			for (auto i2=0; i2 < this->_box_length; i2++) {
				for (auto i3=0; i3 < this->_box_length; i3++) {
					spair = get_mol_it(Site(idxs[0][i1],idxs[1][i2],idxs[2][i3]),sp);
					if (spair.first) 
					{
						return spair;
					};
				};
			};
		};
		return std::make_pair(false,SiteIt());
	};

	/********************
	Get a free site
	********************/

	std::pair<bool,Site> Lattice::get_free_site() 
	{
		// Get random indexes to search
		std::map<int,std::vector<int>> idxs = _get_random_idxs();

		// Try all sites
		Site s;
		std::pair<bool,SiteIt> spair;
		for (auto i1=0; i1 < this->_box_length; i1++) {
			for (auto i2=0; i2 < this->_box_length; i2++) {
				for (auto i3=0; i3 < this->_box_length; i3++) {
					s = Site(idxs[0][i1],idxs[1][i2],idxs[2][i3]);
					spair = get_mol_it(s);
					if (!(spair.first)) {
						return std::make_pair(true,s);
					};
				};
			};
		};

		return std::make_pair(false,s);
	};

	/********************
	Get neighbors of a site
	********************/

	std::pair<Site,std::pair<bool,SiteIt>> Lattice::get_neighbor_random(Site s)
	{
		// All neighbors
		std::vector<Site> nbrs = _get_all_neighbors(s);

		// Shuffle
		std::random_shuffle(nbrs.begin(),nbrs.end());

		// Random choice
		auto it = nbrs.begin();
		std::advance(it,randI(0,nbrs.size()-1));

		// Check if occupied
		std::pair<bool,SiteIt> occ = get_mol_it(*it);
		if (occ.first) {
			// Yes
			return std::make_pair(Site(occ.second),std::make_pair(true,occ.second));
		} else {
			// No
			return std::make_pair(*it,std::make_pair(false,SiteIt()));
		};
	};

	std::pair<Site,std::pair<bool,SiteIt>> Lattice::get_neighbor_random(SiteIt sit)
	{
		return get_neighbor_random(Site(sit));
	};

	std::pair<bool,Site> Lattice::get_free_neighbor_random(Site s) 
	{
		// All allowed nbrs
		std::vector<Site> nbrs = _get_all_neighbors(s);

		// Shuffle
		std::random_shuffle(nbrs.begin(),nbrs.end());

		// Check all neighbors
		std::pair<bool,SiteIt> spair;
		for (auto n: nbrs) {
			spair = get_mol_it(n);
			if (!(spair.first)) {
				return std::make_pair(true,n);
			};
		};

		return std::make_pair(false,Site());
	};

	std::pair<bool,Site> Lattice::get_free_neighbor_random(SiteIt sit) 
	{
		return get_free_neighbor_random(Site(sit));
	};

	/********************
	Get count, NN of species
	********************/

	int Lattice::get_count(std::string sp)
	{
		int c = 0;
		// Go through all sites
		Site s;
		std::pair<bool,SiteIt> spair;
		for (auto i1=1; i1 <= this->_box_length; i1++) {
			for (auto i2=1; i2 <= this->_box_length; i2++) {
				for (auto i3=1; i3 <= this->_box_length; i3++) {
					s = Site(i1,i2,i3);
					spair = get_mol_it(s,sp);
					if (spair.first) {
						c++;
					};
				};
			};
		};
		return c;
	};

	int Lattice::get_nn(std::string sa, std::string sb)
	{
		int nn = 0;
		// Go through all sites
		Site s;
		std::pair<bool,SiteIt> spair;
		std::vector<Site> nbrs;
		for (auto i1=1; i1 <= this->_box_length; i1++) {
			for (auto i2=1; i2 <= this->_box_length; i2++) {
				for (auto i3=1; i3 <= this->_box_length; i3++) {
					s = Site(i1,i2,i3);
					spair = get_mol_it(s,sa);
					if (spair.first) {
						// This one is species A
						// Now search all neigbors
						nbrs = _get_all_neighbors(s);
						// Go through all neighbors
						for (auto nbr: nbrs) {
							spair = get_mol_it(nbr,sb);
							if (spair.first) {
								// Ok!
								nn++;
							};
						};
					};
				};
			};
		};
		// Double counting?
		if (sa == sb) {
			nn /= 2;
		};
		return nn;
	};

	/********************
	Write lattice to a file
	********************/

	void Lattice::write_to_file(std::string fname) 
	{
		std::ofstream f;
		f.open (fname);
		for (auto it = this->_map.begin(); it != this->_map.end(); it++) {
			for (auto it_1 = it->second.begin(); it_1 != it->second.end(); it_1++) {
				for (auto it_2 = it_1->second.begin(); it_2 != it_1->second.end(); it_2++) {
					f << it->first << " " << it_1->first << " " << it_2->first << " " << it_2->second << "\n";
				};
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
			    	make_mol(Site(atoi(x.c_str()),atoi(y.c_str()),atoi(z.c_str())),sp);
			    	//std::cout << "Got: " << atoi(x.c_str()) << " " << atoi(y.c_str()) << " " << atoi(z.c_str()) << " " << sp << std::endl;
			    	i_frag = 0;
			    	sp=""; x=""; y=""; z="";
			    };
			};
		};
		f.close();
	};

	/********************
	Anneal
	********************/

	void Lattice::anneal(std::map<std::string,double> &h_dict,std::map<StrPair,double> &j_dict, int n_steps) {

		Site s;

		std::pair<bool,SiteIt> ret_site;
		std::string sp;

		std::vector<Site> nbrs;
		std::vector<std::string> nbrs_occ;
		std::vector<Site>::iterator it_nbr;
		std::pair<bool,SiteIt> ret_nbr;

		std::map<std::string,double>::iterator ith;

		double hOld,jOld,hNew,jNew,energy_diff;

		// Go through the steps
		for (int i=0; i<n_steps; i++) {

			// Pick a site to flip randomly
			s = Site(randI(1,_box_length),randI(1,_box_length),randI(1,_box_length));

			// Get occupied neighbors 
			nbrs = _get_all_neighbors(s);
			it_nbr = nbrs.begin();
			while (it_nbr != nbrs.end()) {
				ret_nbr = get_mol_it(*it_nbr);
				if (ret_nbr.first) {
					// Occupied
					nbrs_occ.push_back(ret_nbr.second.it_2->second);
					it_nbr++;
				} else {
					// Empty
					it_nbr = nbrs.erase(it_nbr);
				};
			};

			// Check if this site is occupied
			ret_site = get_mol_it(s);
			if (ret_site.first) {
				// Occupied, flip down
				sp = ret_site.second.it_2->second;
				hOld = -h_dict[sp];
				jOld = 0.0;
				for (auto nbr: nbrs_occ) {
					jOld -= j_dict[StrPair(sp,nbr)];
				};
				// New couplings
				hNew = 0.0;
				jNew = 0.0;
			} else {
				// Unoccupied, flip up
				hOld = 0.0;
				jOld = 0.0;
				// Random species
				ith = h_dict.begin();
				std::advance(ith, randI(0,h_dict.size()-1));
				sp = ith->first;
				// New couplings
				hNew = -ith->second;
				jNew = 0.0;
				for (auto nbr: nbrs_occ) {
					jNew -= j_dict[StrPair(sp,nbr)];
				};
			};

			// Energy difference
			energy_diff = hNew + jNew - hOld - jOld;
			if (energy_diff < 0.0 || exp(-energy_diff) > randD(0.0,1.0)) {
				// Accept the flip!
				if (ret_site.first) {
					// Occupied, flip down
					erase_mol_it(ret_site.second);
				} else {
					// Unoccupied, flip up
					make_mol(s,sp);
				};
			};
		};
	};

	/****************************************
	Lattice - PRIVATE
	****************************************/

	/********************
	Get random indexes
	********************/

	std::map<int,std::vector<int>> Lattice::_get_random_idxs()
	{
		// Vectors of sites
		std::vector<int> x1(this->_box_length);
		std::iota(std::begin(x1), std::end(x1), 1);
		std::vector<int> x2 = x1, x3 = x1;

		// Shuffle
		std::random_shuffle(x1.begin(),x1.end());
		std::random_shuffle(x2.begin(),x2.end());
		std::random_shuffle(x3.begin(),x3.end());

		// Return
		std::map<int,std::vector<int>> m;
		m[0] = x1;
		m[1] = x2;
		m[2] = x3;
		return m;
	};

	/********************
	Get all neighbors of a site
	********************/

	std::vector<Site> Lattice::_get_all_neighbors(Site s)
	{
		std::vector<Site> nbrs;

		// 6 are possible
		Site nbr_Site = s;

		// x
		nbr_Site.x += 1;
		if (nbr_Site.x <= this->_box_length) {
			nbrs.push_back(nbr_Site);
		};
		nbr_Site.x -= 2;
		if (nbr_Site.x >= 1) {
			nbrs.push_back(nbr_Site);
		};
		nbr_Site.x += 1;
		// y
		nbr_Site.y += 1;
		if (nbr_Site.y <= this->_box_length) {
			nbrs.push_back(nbr_Site);
		};
		nbr_Site.y -= 2;
		if (nbr_Site.y >= 1) {
			nbrs.push_back(nbr_Site);
		};
		nbr_Site.y += 1;
		// z
		nbr_Site.z += 1;
		if (nbr_Site.z <= this->_box_length) {
			nbrs.push_back(nbr_Site);
		};
		nbr_Site.z -= 2;
		if (nbr_Site.z >= 1) {
			nbrs.push_back(nbr_Site);
		};
		nbr_Site.z += 1;

		return nbrs;
	};

};