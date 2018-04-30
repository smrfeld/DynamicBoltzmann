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
	Structure to hold pairs and triplets of sites
	****************************************/

	class Site;

	struct Site2 {
		Site *s1,*s2;
		Site2(Site* s1, Site *s2) {
			this->s1 = s1;
			this->s2 = s2;
		};
	};

	struct Site3 {
		Site *s1,*s2,*s3;
		Site3(Site* s1, Site *s2, Site *s3) {
			this->s1 = s1;
			this->s2 = s2;
			this->s3 = s3;
		};
	};

	/****************************************
	Class to hold a connection from visible to hidden
	****************************************/

	class HiddenUnit;

	class ConnectionVH {
	private:

		// Site
		Site *_site;

		// Hidden unit
		HiddenUnit *_hidden_unit;

		// Ixn param W associated with this connection
		// Stored both directions
		std::map<Species*, std::map<HiddenSpecies*, std::vector<IxnParam*> > > _ips_visible_hidden;
		std::map<HiddenSpecies*, std::map<Species*, std::vector<IxnParam*> > > _ips_hidden_visible;

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const ConnectionVH& other);

	public:

		/********************
		Constructor
		********************/

		ConnectionVH(Site *site, HiddenUnit *hidden_unit, std::vector<IxnParam*> ips);
		ConnectionVH(const ConnectionVH& other);
		ConnectionVH(ConnectionVH&& other);
		ConnectionVH& operator=(const ConnectionVH& other);
		ConnectionVH& operator=(ConnectionVH&& other);
		~ConnectionVH();

		/********************
		Add ixn param
		********************/

		void add_ixn_param(IxnParam* ip);

		/********************
		Get activation for a species on the visible unit
		********************/

		double get_act_visible(Species* sp_visible);

		/********************
		Get activation for a species on the hidden unit
		********************/

		double get_act_hidden(HiddenSpecies* sp_hidden);
	};

	/****************************************
	Class to hold a lattice site
	***************************************/

	class Site {
	private:

		// Dimensionality and location
		int _dim;
		int _x;
		int _y;
		int _z;

		// Probs of species/empty
		std::vector<Species*> _sp_possible;
		std::map<Species*, double> _probs;
		double _prob_empty;

		// Connectivity to any hidden units
		// A species-dependent graph :)
		std::vector<ConnectionVH*> _hidden_conns;

		// Neighbors
		std::vector<Site*> _nbrs;
		std::vector<Site2> _nbrs_triplets;
		std::vector<Site3> _nbrs_quartics;

		// Add/Remove counts on a given species (not _sp_binary, unless it is passed)
		void _remove_counts_on_species(Species *sp, double prob);
		void _add_counts_on_species(Species *sp, double prob);

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const Site& other);

	public:

		/********************
		Constructor
		********************/

		Site(int x);
		Site(int x, int y);
		Site(int x, int y, int z);
		Site(const Site& other);
		Site(Site&& other);
		Site& operator=(const Site& other);
		Site& operator=(Site&& other);
		~Site();

		/********************
		Check location
		********************/

		int x() const;
		int y() const;
		int z() const;
		bool less_than(const Site &other) const;

		/********************
		Add neighbors
		********************/

		void add_nbr(Site *s);
		void add_nbr_triplet(Site *s1, Site *s2);
		void add_nbr_quartic(Site *s1, Site *s2, Site *s3);

		/********************
		Add a hidden conn
		********************/

		void add_visible_hidden_conn(ConnectionVH* connvh);

		/********************
		Add a species possibility
		********************/

		void add_species_possibility(Species* sp);

		/********************
		Get probability
		********************/

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

		/********************
		Sample
		********************/

		void sample(bool binary);

		/********************
		Check if site is empty
		********************/

		bool empty() const;

		/********************
		Binarize the set
		********************/

		void binarize();
	};
	// Comparator
	bool operator <(const Site& a, const Site& b);

	/****************************************
	Lattice
	****************************************/

	typedef std::list<Site> lattice;
	typedef std::list<Site>::iterator latt_it;

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
		Site* _look_up(int x);
		Site* _look_up(int x, int y);
		Site* _look_up(int x, int y, int z);

		// Pointers to species present
		std::map<std::string,Species*> _sp_map;
		std::vector<Species*> _sp_vec;

		// Does the lattice have the following structure... ?
		bool _latt_has_nn_structure;
		bool _latt_has_triplet_structure;
		bool _latt_has_quartic_structure;

		// Contructor helpers
		void _clean_up();
		void _reset();
		void _copy(const Lattice& other);

	public:

		/********************
		Constructor
		********************/

		Lattice(int dim, int box_length);
		Lattice(const Lattice& other);
		Lattice(Lattice&& other);
		Lattice& operator=(const Lattice& other);
		Lattice& operator=(Lattice&& other);
		~Lattice();

		/********************
		Getters
		********************/

		int dim() const;
		int box_length() const;

		/********************
		Add a species
		********************/

		void add_species_possibility(Species *sp);

		/********************
		Initialize lattice structure of NNs, triplets, etc
		********************/

		void init_nn_structure();
		void init_triplet_structure();
		void init_quartic_structure();

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
		Binarize
		********************/

		void binarize();

		/********************
		Write lattice to a file
		********************/

		void write_to_file(std::string fname);

		/********************
		Read lattice from a file
		********************/

		void read_from_file(std::string fname, bool binary=true);

		/********************
		Populate randomly according to some counts
		********************/

		void populate_randomly();
		void populate_randomly(std::map<Species*, int> counts);

		/********************
		Sample
		********************/

		void sample(bool binary=true);

	};

};