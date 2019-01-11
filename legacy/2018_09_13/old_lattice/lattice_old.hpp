// string
#ifndef STRING_h
#define STRING_h
#include <string>
#endif

// vector
#ifndef VECTOR_h
#define VECTOR_h
#include <vector>
#endif

// map
#ifndef MAP_h
#define MAP_h
#include <map>
#endif

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	General functions
	****************************************/

	typedef std::map<int,std::map<int,std::map<int,std::string>>> lattice_map;
	typedef std::map<int,std::map<int,std::string>> lattice_map_1;
	typedef std::map<int,std::string> lattice_map_2;

	/****************************************
	Structure to hold a lattice site iterator
	****************************************/

	struct SiteIt {
		lattice_map::iterator it;
		lattice_map_1::iterator it_1;
		lattice_map_2::iterator it_2;	

		// Constructor
		SiteIt();
		SiteIt(lattice_map::iterator itIn, lattice_map_1::iterator it_1In, lattice_map_2::iterator it_2In);	
	};
	std::ostream& operator<<(std::ostream& os, const SiteIt& sit);

	/****************************************
	Structure to hold a lattice site
	****************************************/

	struct Site {
		int x;
		int y;
		int z;	

		// Constructor
		Site();
		Site(int xIn, int yIn, int zIn);
		Site(SiteIt sit);
	};
	// Comparator
	bool operator <(const Site& a, const Site& b);
	bool operator==(const Site& a, const Site& b);
	std::ostream& operator<<(std::ostream& os, const Site& s);

	/****************************************
	Lattice
	****************************************/

	class Lattice
	{
	private:

		// Internal maps
		lattice_map _map;

		// Size
		int _box_length;

		/********************
		Get random indexes
		********************/

		std::map<int,std::vector<int>> _get_random_idxs();

		/********************
		Get all neighbors of a site
		********************/

		std::vector<Site> _get_all_neighbors(Site s);

	public:

		/********************
		Constructor/Destructor
		********************/

		Lattice(int box_length);
		Lattice();
		~Lattice();

		/********************
		Clear, size
		********************/

		void clear();
		int size();

		/********************
		Make a mol
		********************/

		std::pair<bool,SiteIt> make_mol(Site s, std::string sp);
		std::pair<bool,SiteIt> make_mol_random(std::string sp);

		/********************
		Erase a mol
		********************/

		bool erase_mol(Site s);
		bool erase_mol_it(SiteIt sit);
		std::pair<bool,Site> erase_mol_random(std::string sp);

		/********************
		Get a mol
		********************/

		std::pair<bool,SiteIt> get_mol_it(Site s);
		std::pair<bool,SiteIt> get_mol_it(Site s, std::string sp);
		std::pair<bool,SiteIt> get_mol_random_it();
		std::pair<bool,SiteIt> get_mol_random_it(std::string sp);

		/********************
		Get a free site
		********************/

		std::pair<bool,Site> get_free_site();

		/********************
		Get neighbors of a site
		********************/

		std::pair<Site,std::pair<bool,SiteIt>> get_neighbor_random(Site s);
		std::pair<Site,std::pair<bool,SiteIt>> get_neighbor_random(SiteIt sit);
		std::pair<bool,Site> get_free_neighbor_random(Site s);
		std::pair<bool,Site> get_free_neighbor_random(SiteIt sit);

		/********************
		Get count,NN of species
		********************/

		int get_count(std::string sp);
		int get_nn(std::string sa, std::string sb);

		/********************
		Write lattice to a file
		********************/

		void write_to_file(std::string fname);

		/********************
		Read lattice from a file
		********************/

		void read_from_file(std::string fname);

		/********************
		Anneal
		********************/

		void anneal(std::map<std::string,double> &h_dict,std::map<std::string,std::map<std::string,double>> &j_dict, int n_steps);
	};

};