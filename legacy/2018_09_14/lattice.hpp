#ifndef STRING_H
#define STRING_H
#include <string>
#endif

#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

#ifndef LIST_H
#define LIST_H
#include <list>
#endif

#ifndef MAP_H
#define MAP_H
#include <map>
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	// Forwards
	class Site;
	class Species;
	struct Species2;
	struct Species3;
	class IxnParamTraj;

	/****************************************
	Lattice
	****************************************/

	typedef std::list<Site> Latt;
	typedef std::list<Site>::iterator Latt_it;

	class Lattice
	{
	private:

		// Dimensionality
		int _dim;

		// Size
		int _box_length;

		// Internal maps
		Latt _latt;

		// Lookup a site iterator from x,y,z
		Site* _look_up(int x);
		Site* _look_up(int x, int y);
		Site* _look_up(int x, int y, int z);

		// Pointers to species present
		std::map<std::string,Species*> _sp_map;
		std::vector<Species*> _sp_vec;

		// Does the Latt have the following structure... ?
		bool _latt_has_nn_structure;
		bool _latt_has_triplet_structure;
		bool _latt_has_quartic_structure;

		// Counts
		std::map<Species*, int> _counts;
		std::map<Species2, int> _nn_counts;
		std::map<Species3, int> _triplet_counts;
		std::map<Species4, int> _quartic_counts;

		// Ixns
		std::map<Species*, std::vector<IxnParamTraj*>> _h_ptrs;
		std::map<Species2, std::vector<IxnParamTraj*>> _j_ptrs;
		std::map<Species3, std::vector<IxnParamTraj*>> _k_ptrs;

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
		Validate graph
		********************/

		void validate_graph() const;

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
		Initialize Latt structure of NNs, triplets, etc
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
		Write Latt to a file
		********************/

		void write_to_file(std::string fname);

		/********************
		Read Latt from a file
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