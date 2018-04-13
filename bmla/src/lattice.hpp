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

		// Probs of species/empty
		std::map<Species*, double> _probs;
		double _prob_empty;

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const Site& other);

		// Add/Remove counts on a given species (not _sp_binary, unless it is passed)
		void _remove_counts_on_species(Species *sp, double prob);
		void _add_counts_on_species(Species *sp, double prob);

	public:
		int dim;
		int x;
		int y;
		int z;
		std::vector<latt_it> nbrs;
		std::vector<std::pair<latt_it,latt_it>> triplets;

		// Connectivity to any hidden units
		// A species-dependent graph :)
		std::map<Species*, std::vector<HiddenUnit*>> hidden_conns;

		// Constructor
		Site();
		Site(int xIn);
		Site(int xIn, int yIn);
		Site(int xIn, int yIn, int zIn);
		Site(const Site& other);
		Site(Site&& other);
		Site& operator=(const Site& other);
		Site& operator=(Site&& other);
		~Site();

		// Add a species possibility
		void add_species(Species* sp);

		// Get a probability
		// If binary, returns 1.0 for a certain species
		// Pass nullptr to get prob of empty
		double get_prob(Species *sp) const;
		// Get all probs
		const std::map<Species*, double>& get_probs() const;
		// Set probability
		// Pass nullptr to set probability of being empty
		void set_prob(Species *sp, double prob);
		// Set site to be empty
		void set_site_empty();
		// Set site to have binary probs
		void set_site_binary(Species *sp);

		// Get J and K activations
		// Go through all possible species: (probabilities) times (J of the coupling for the given species)
		double get_act_j(Species *sp) const;
		double get_act_k(Species *sp) const;

		// Check if site is empty
		bool empty() const;
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
		bool _exists_hidden;
		bool _exists_h;
		bool _exists_j;
		bool _exists_k;

		// Sample an unnormalized probability vector
		int _sample_prop_vec(std::vector<double> &props);
		
		// Sample a vector of propensities (cumulative probabilities)
		int _sample_prob_vec(std::vector<double> &probs);

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

		void set_exists_hidden(bool flag=true);
		void set_exists_h(bool flag=true);
		void set_exists_j(bool flag=true);
		void set_exists_k(bool flag=true);

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

		void sample(bool binary=true);

	};

};