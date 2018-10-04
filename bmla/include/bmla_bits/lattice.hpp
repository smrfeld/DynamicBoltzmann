#include <string>
#include <vector>
#include <map>

#ifndef FWDS_SPECIES_H
#define FWDS_SPECIES_H
#include "fwds/fwds_species.hpp"
#endif

/************************************
* Namespace for bmla
************************************/

namespace bmla {

	// Forwards
	class UnitVisible;
	class UnitHidden;
	class BiasDict;
	class ConnVH;
	class ConnVV;
	class ConnVVV;
	class O2IxnDict;
	class O3IxnDict;
	class Moment;

	/****************************************
	Lattice
	****************************************/

	class Lattice
	{
	private:

		// Dimensionality
		int _no_dims;

		// Size
		int _box_length;

		// Visible layer
		std::vector<UnitVisible*> _latt_v;

		// Hidden layers
		std::map<int,std::vector<UnitHidden*>> _latt_h;

		// Connections
		std::vector<ConnVV*> _conns_vv;
		std::vector<ConnVVV*> _conns_vvv;
		std::vector<ConnVH*> _conns_vh;

		// Lookup a site iterator from x,y,z
		UnitVisible* _look_up_unit_v(int x) const;
		UnitVisible* _look_up_unit_v(int x, int y) const;
		UnitVisible* _look_up_unit_v(int x, int y, int z) const;
		UnitHidden* _look_up_unit_h(int layer, int idx) const;

		// Check dim
		void _check_dim(int dim) const;

		// Count helpers
		void _get_count(double &count, Sptr &sp1, Sptr &sp2, const UnitVisible *uv1, const UnitVisible *uv2, bool binary, bool reversibly) const;
		void _get_count(double &count, Sptr &sp1, Sptr &sp2, Sptr &sp3, const UnitVisible *uv1, const UnitVisible *uv2, const UnitVisible *uv3, bool binary, bool reversibly) const;

		// Contructor helpers
		void _clean_up();
		void _move(Lattice& other);
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
		Check setup
		********************/

		void print() const;

		/********************
		Helpers to setup all sites
		********************/

		// Add possible species
		void all_units_v_add_possible_species(Sptr species);
		void all_units_h_add_possible_species(Sptr species);

		// Biases
		void all_units_v_set_bias_dict(std::shared_ptr<BiasDict> bias_dict);
		void all_units_h_set_bias_dict(std::shared_ptr<BiasDict> bias_dict);

		// Make connections
		void all_conns_vv_init();
		void all_conns_vv_init(std::shared_ptr<O2IxnDict> ixn_dict);
		void all_conns_vvv_init();
		void all_conns_vvv_init(std::shared_ptr<O3IxnDict> ixn_dict);

		// Set ixn dicts of connections
		void all_conns_vv_set_ixn_dict(std::shared_ptr<O2IxnDict> ixn_dict);
		void all_conns_vvv_set_ixn_dict(std::shared_ptr<O3IxnDict> ixn_dict);
		void all_conns_vh_set_ixn_dict(std::shared_ptr<O2IxnDict> ixn_dict);

		// Link units to moments
		void all_units_v_add_to_moment_h(std::shared_ptr<Moment> moment);
		void all_units_h_add_to_moment_b(std::shared_ptr<Moment> moment);
		void all_conns_vv_add_to_moment_j(std::shared_ptr<Moment> moment);
		void all_conns_vvv_add_to_moment_k(std::shared_ptr<Moment> moment);
		void all_conns_vh_add_to_moment_w(std::shared_ptr<Moment> moment);

		/********************
		Add visible-visible connections
		********************/

		void add_conn_vv(UnitVisible *uv1, UnitVisible *uv2);
		void add_conn_vv(UnitVisible *uv1, UnitVisible *uv2, std::shared_ptr<O2IxnDict> ixn_dict);
		void add_conn_vvv(UnitVisible *uv1, UnitVisible *uv2, UnitVisible *uv3);
		void add_conn_vvv(UnitVisible *uv1, UnitVisible *uv2, UnitVisible *uv3, std::shared_ptr<O3IxnDict> ixn_dict);

		/********************
		Add hidden units
		********************/

		// Returns the idx of this unit
		int add_hidden_unit(int layer);
		int add_hidden_unit(int layer, std::vector<Sptr> species_possible);

		/********************
		Add visible-hidden connections
		********************/

		void add_conn_vh(UnitVisible *uv, UnitHidden *uh);
		void add_conn_vh(UnitVisible *uv, UnitHidden *uh, std::shared_ptr<O2IxnDict> ixn_dict);

		/********************
		Get unit
		********************/

		UnitVisible* get_unit_v(int x) const;
		UnitVisible* get_unit_v(int x, int y) const;
		UnitVisible* get_unit_v(int x, int y, int z) const;

		UnitHidden* get_unit_h(int layer, int idx) const;

		/********************
		Get connection
		********************/

		ConnVV* get_conn_vv(int x1, int x2) const;
		ConnVV* get_conn_vv(int x1, int y1, int x2, int y2) const;
		ConnVV* get_conn_vv(int x1, int y1, int z1, int x2, int y2, int z2) const;

		ConnVVV* get_conn_vvv(int x1, int x2, int x3) const;
		ConnVVV* get_conn_vvv(int x1, int y1, int x2, int y2, int x3, int y3) const;
		ConnVVV* get_conn_vvv(int x1, int y1, int z1, int x2, int y2, int z2, int x3, int y3, int z3) const;

		ConnVH* get_conn_vh(int x, int layer, int idx) const;
		ConnVH* get_conn_vh(int x, int y, int layer, int idx) const;
		ConnVH* get_conn_vh(int x, int y, int z, int layer, int idx) const;

		/********************
		Getters (general)
		********************/

		int get_no_dims() const;
		int get_no_units_v();
		int get_no_units_h();

		/********************
		Apply funcs to all units
		********************/

		// Clear the lattice
		void all_units_v_set_empty();
		void all_units_h_set_empty();
		void all_units_set_empty();

		// Binary/probabilistic
		void all_units_convert_to_b_mode();
		void all_units_convert_to_p_mode();

		/********************
		Write/read latt to a file
		********************/

		void write_to_file(std::string fname, bool binary=true);
		void read_from_file(std::string fname, bool binary=true);

		/********************
		Sample
		********************/

		void sample_v(bool binary=true);
		void sample_h(bool binary=true);

		/********************
		Get counts
		********************/

		double get_count(Sptr &sp, bool binary=true) const;
		double get_count(Sptr &sp1, Sptr &sp2, bool binary=true, bool reversibly=true) const;
		double get_count(Sptr &sp1, Sptr &sp2, Sptr &sp3, bool binary=true, bool reversibly=true) const;
		double get_count(Sptr &sp1, Sptr &sp2, Sptr &sp3, Sptr &sp4, bool binary=true, bool reversibly=true) const;
	};

};