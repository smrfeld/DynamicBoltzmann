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
	class ConnHH;
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
		std::vector<int> _latt_v_idxs;

		// Hidden layers
		std::map<int,std::vector<UnitHidden*>> _latt_h;
		std::map<int,std::vector<int>> _latt_h_idxs;

		// Connections
		std::vector<ConnHH*> _conns_hh;
		std::vector<ConnVV*> _conns_vv;
		std::vector<ConnVVV*> _conns_vvv;
		std::vector<ConnVH*> _conns_vh;

		// IO
		std::map<std::string,Sptr> _IO_species_possible;
		bool _IO_did_init;

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

		ConnVV* add_conn_vv(UnitVisible *uv1, UnitVisible *uv2);
		ConnVV* add_conn_vv(UnitVisible *uv1, UnitVisible *uv2, std::shared_ptr<O2IxnDict> ixn_dict);
		ConnVVV* add_conn_vvv(UnitVisible *uv1, UnitVisible *uv2, UnitVisible *uv3);
		ConnVVV* add_conn_vvv(UnitVisible *uv1, UnitVisible *uv2, UnitVisible *uv3, std::shared_ptr<O3IxnDict> ixn_dict);

		/********************
		Add hidden units
		********************/

		UnitHidden* add_hidden_unit(int layer);
		UnitHidden* add_hidden_unit(int layer, std::vector<Sptr> species_possible);

		/********************
		Add visible-hidden connections
		********************/

		ConnVH* add_conn_vh(UnitVisible *uv, UnitHidden *uh);
		ConnVH* add_conn_vh(UnitVisible *uv, UnitHidden *uh, std::shared_ptr<O2IxnDict> ixn_dict);

		/********************
		Add hidden-hidden connections
		********************/

		ConnHH* add_conn_hh(UnitHidden *uh1, int layer_1, UnitHidden *uh2, int layer_2);
		ConnHH* add_conn_hh(UnitHidden *uh1, int layer_1, UnitHidden *uh2, int layer_2, std::shared_ptr<O2IxnDict> ixn_dict);

		/********************
		Get unit
		********************/

		const std::vector<UnitVisible*>& get_all_units_v() const;
		const std::vector<UnitHidden*>& get_all_units_h(int layer) const;
		const std::map<int,std::vector<UnitHidden*>>& get_all_units_h() const;

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

		ConnHH* get_conn_hh(int layer1, int idx1, int layer2, int idx2) const;

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

		// Random
		void all_units_v_random();

		// Binarize
		void all_units_v_binarize();
		void all_units_h_binarize();

		/********************
		Write/read latt to a file
		********************/

		void write_to_file(std::string fname, bool binary);

		void init_file_reader(std::vector<Sptr> species_possible);
		void read_from_file(std::string fname, bool binary);

		/********************
		Sample
		********************/

		void sample_down_h_to_v(bool layer_wise, bool binary_visible, bool binary_hidden);
		void sample_up_v_to_h(bool layer_wise, bool binary_hidden);

		/********************
		Get counts
		********************/

		double get_count(Sptr &sp) const;
		double get_count(Sptr &sp1, Sptr &sp2, bool reversibly) const;
		double get_count(Sptr &sp1, Sptr &sp2, Sptr &sp3, bool reversibly) const;
		double get_count(Sptr &sp1, Sptr &sp2, Sptr &sp3, Sptr &sp4, bool reversibly) const;
	};

};