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

// list
#ifndef LIST_h
#define LIST_h
#include <list>
#endif

// map
#ifndef MAP_h
#define MAP_h
#include <map>
#endif

// Species
#ifndef SPECIES_h
#define SPECIES_h
#include "species.hpp"
#endif

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	General functions
	****************************************/

	class Site;
	typedef std::list<Site> lattice;
	typedef std::list<Site>::iterator latt_it;

	/****************************************
	Structure to hold a lattice site
	****************************************/

	class Site {
	private:

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const Site& other);

	public:
		int dim;
		int x;
		int y;
		int z;	
		Species *sp;
		std::vector<latt_it> nbrs;

		// Connectivity to any hidden units
		// A species-dependent graph :)
		std::map<Species*, std::vector<HiddenUnit*>> hidden_conns;

		// Constructor
		Site();
		Site(int xIn);
		Site(int xIn, Species *spIn);
		Site(int xIn, int yIn);
		Site(int xIn, int yIn, Species *spIn);
		Site(int xIn, int yIn, int zIn);
		Site(int xIn, int yIn, int zIn, Species *spIn);
		Site(const Site& other);
		Site(Site&& other);
		Site& operator=(const Site& other);
		Site& operator=(Site&& other);
		~Site();
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

		// Dimensionality
		int _dim;

		// Size
		int _box_length;

		// Internal maps
		lattice _latt;

		// Lookup a site iterator from x,y,z
		latt_it _look_up(int x);
		latt_it _look_up(int x, int y);
		latt_it _look_up(int x, int y, int z);

		// Pointers to species present
		std::map<std::string,Species*> _sp_map;
		std::vector<Species*> _sp_vec;

		// Flag - are hidden units present? (needed for annealing)
		bool _hidden_layer_exists;

		// Contructor helpers
		void _clean_up();
		void _reset();
		void _copy(const Lattice& other);

	public:

		/********************
		Constructor
		********************/

		Lattice(int dim, int box_length);
		Lattice();
		Lattice(const Lattice& other);
		Lattice(Lattice&& other);
		Lattice& operator=(const Lattice& other);
		Lattice& operator=(Lattice&& other);
		~Lattice();

		/********************
		Getters
		********************/

		int dim() const;

		/********************
		Add a species
		********************/

		void add_species(Species *sp);

		/********************
		Indicate that the hidden unit exists
		********************/

		void set_hidden_layer_exists();

		/********************
		Find a pointer to a site by index
		********************/

		Site* get_site(int x);
		Site* get_site(int x, int y);
		Site* get_site(int x, int y, int z);

		/********************
		Clear, size
		********************/

		void clear();
		int size();

		/********************
		Make a mol
		********************/

		bool make_mol(latt_it s, Species *sp);

		/********************
		Erase a mol
		********************/

		bool erase_mol(latt_it s);

		/********************
		Write lattice to a file
		********************/

		void write_to_file(std::string fname);

		/********************
		Read lattice from a file
		********************/

		void read_from_file(std::string fname);

		/********************
		Sample
		********************/

		void sample();

		/********************
		Sample probabilities/propensities
		********************/

		// Sample an unnormalized probability vector
		int sample_prop_vec(std::vector<double> &props);
		
		// Sample a vector of propensities (cumulative probabilities)
		int sample_prob_vec(std::vector<double> &probs);
	};

};