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

#ifndef UNORDERED_MAP_H
#define UNORDERED_MAP_H
#include <unordered_map>
#endif

#ifndef FWDS_SPECIES_H
#define FWDS_SPECIES_H
#include "fwds/fwds_species.hpp"
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	// Forwards
	class UnitVisible;
	class BiasDict;
	class ConnVV;
	class ConnVVV;
	class O2IxnDict;
	class O3IxnDict;
	struct Sptr2;
	struct Sptr3;

	/****************************************
	Lattice
	****************************************/

	typedef std::list<UnitVisible> Latt;
	typedef Latt::iterator Latt_it;

	class Lattice
	{
	private:

		// Dimensionality
		int _dim;

		// Size
		int _box_length;

		// Visible units
		Latt _latt;

		// Connections
		std::list<ConnVV> _conns_vv;
		std::list<ConnVVV> _conns_vvv;

		// Lookup a site iterator from x,y,z
		UnitVisible& _look_up(int x);
		UnitVisible& _look_up(int x, int y);
		UnitVisible& _look_up(int x, int y, int z);

		// Check dim
		void _check_dim(int dim) const;

		// Contructor helpers
		void _clean_up();
		void _reset();
		void _copy(const Lattice& other);

	public:

		/********************
		Constructor
		********************/

		Lattice(int dim, int box_length);
		Lattice(int dim, int box_length, std::vector<Sptr> possible_species_of_all_units_visible);
		Lattice(const Lattice& other);
		Lattice(Lattice&& other);
		Lattice& operator=(const Lattice& other);
		Lattice& operator=(Lattice&& other);
		~Lattice();

		/********************
		Check setup
		********************/

		void check_setup() const;
		void print_occupancy(bool binary=true) const;

		/********************
		Finish setup
		********************/

		// Add possible species
		void add_possible_species_to_all_units_visible(Sptr species);
		// void add_possible_species_to_all_unit_hidden(Sptr species);

		// Biases
		void set_bias_dict_of_all_units_visible(std::shared_ptr<BiasDict> bias_dict);
		// void set_bias_dict_of_all_units_hidden(std::shared_ptr<BiasDict> bias_dict);

		// Visible-Visible ixns

		void init_conns_NN_all_units_visible();
		void init_conns_NN_all_units_visible(std::shared_ptr<O2IxnDict> ixn_dict);

		void init_conns_triplet_all_units_visible();
		void init_conns_triplet_all_units_visible(std::shared_ptr<O3IxnDict> ixn_dict);

		void set_ixn_dict_all_conns_vv(std::shared_ptr<O2IxnDict> ixn_dict);
		void set_ixn_dict_all_conns_vvv(std::shared_ptr<O3IxnDict> ixn_dict);

		void add_conn_vv(UnitVisible *uv1, UnitVisible *uv2, std::shared_ptr<O2IxnDict> ixn_dict);
		void add_conn_vvv(UnitVisible *uv1, UnitVisible *uv2, UnitVisible *uv3, std::shared_ptr<O3IxnDict> ixn_dict);

		// Hidden units
		// void add_hidden_unit();

		/********************
		Get unit
		********************/

		UnitVisible& get_unit_visible(int x);
		UnitVisible& get_unit_visible(int x, int y);
		UnitVisible& get_unit_visible(int x, int y, int z);

		/*
		UnitHidden& get_unit_hidden(int layer, int x) const;
		UnitHidden& get_unit_hidden(int layer, int x, int y) const;
		UnitHidden& get_unit_hidden(int layer, int x, int y, int z) const;
		*/

		/********************
		Get connection
		********************/

		ConnVV& get_conn_vv(int x1, int x2);
		ConnVV& get_conn_vv(int x1, int y1, int x2, int y2);
		ConnVV& get_conn_vv(int x1, int y1, int z1, int x2, int y2, int z2);

		ConnVVV& get_conn_vvv(int x1, int x2, int x3);
		ConnVVV& get_conn_vvv(int x1, int y1, int x2, int y2, int x3, int y3);
		ConnVVV& get_conn_vvv(int x1, int y1, int z1, int x2, int y2, int z2, int x3, int y3, int z3);

		/********************
		Getters (general)
		********************/

		int dim() const;
		int no_units_visible();
		// int no_units_hidden();

		/********************
		Apply funcs to all units
		********************/

		// Clear the lattice
		void set_all_units_empty();

		// Binary/probabilistic
		void convert_all_units_to_b_mode();
		void convert_all_units_to_p_mode();

		/********************
		Write/read latt to a file
		********************/

		void write_to_file(std::string fname);
		void read_from_file(std::string fname, bool binary=true);

		/********************
		Sample
		********************/

		void sample_at_timepoint(int timepoint, bool binary=true);

		/********************
		Get counts
		********************/

		double get_count(Sptr &sp, bool binary=true) const;
		double get_count(Sptr &sp1, Sptr &sp2, bool binary=true, bool this_order=false) const;
		double get_count(Sptr &sp1, Sptr &sp2, Sptr &sp3, bool binary=true, bool this_order=false) const;
	};

};