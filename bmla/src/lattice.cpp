#include "../include/bmla_bits/lattice.hpp"

// Other headers
#include "../include/bmla_bits/general.hpp"
#include "../include/bmla_bits/species.hpp"
#include "../include/bmla_bits/unit.hpp"
#include "../include/bmla_bits/connections.hpp"
#include "../include/bmla_bits/ixn_dicts.hpp"
#include "../include/bmla_bits/moment.hpp"

#include <iostream>
#include <fstream>
#include <numeric>
#include "math.h"
#include <ctime>
#include <sstream>
#include <random>

/************************************
* Namespace for bmla
************************************/

namespace bmla {

	/****************************************
	Lattice
	****************************************/

	/********************
	Constructor
	********************/

	Lattice::Lattice(int dim, int box_length)
	{
		if (dim != 1 && dim != 2 && dim != 3) {
			std::cerr << "ERROR: only dimensions 1,2,3 are supported for Lattice." << std::endl;
			exit(EXIT_FAILURE);
		};
		_no_dims = dim;
		_box_length = box_length;

		// Make a fully linked list of sites
		if (dim == 1) {
			for (int x=1; x<=box_length; x++) {
				_latt_v.push_back(new UnitVisible(x));
				_latt_v_idxs.push_back(_latt_v.size()-1);					
			};
		} else if (dim == 2) {
			for (int x=1; x<=box_length; x++) {
				for (int y=1; y<=box_length; y++) {
					_latt_v.push_back(new UnitVisible(x,y));
					_latt_v_idxs.push_back(_latt_v.size()-1);									
				};
			};
		} else if (dim == 3) {
			for (int x=1; x<=box_length; x++) {
				for (int y=1; y<=box_length; y++) {
					for (int z=1; z<=box_length; z++) {
						_latt_v.push_back(new UnitVisible(x,y,z));		
						_latt_v_idxs.push_back(_latt_v.size()-1);								
					};
				};
			};
		};
	};
	Lattice::Lattice(const Lattice& other) {
		_copy(other);
	};
	Lattice::Lattice(Lattice&& other) {
		_move(other);
	};
	Lattice& Lattice::operator=(const Lattice& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	Lattice& Lattice::operator=(Lattice&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
	};
	Lattice::~Lattice() {
		_clean_up();
	};

	void Lattice::_clean_up() {
		for (auto i=0; i<_latt_v.size(); i++) {
			delete _latt_v[i];
			_latt_v[i] = nullptr;
		};
		for (auto &pr: _latt_h) {
			for (auto i=0; i<pr.second.size(); i++) {
				delete pr.second[i];
				pr.second[i] = nullptr;
			};
		};
		for (auto i=0; i<_conns_hh.size(); i++) {
			delete _conns_hh[i];
			_conns_hh[i] = nullptr;
		};
		for (auto i=0; i<_conns_vv.size(); i++) {
			delete _conns_vv[i];
			_conns_vv[i] = nullptr;
		};
		for (auto i=0; i<_conns_vvv.size(); i++) {
			delete _conns_vvv[i];
			_conns_vvv[i] = nullptr;
		};
		for (auto i=0; i<_conns_vh.size(); i++) {
			delete _conns_vh[i];
			_conns_vh[i] = nullptr;
		};
	};
	void Lattice::_copy(const Lattice& other) {
		_no_dims = other._no_dims;
		_box_length = other._box_length;
		for (auto i=0; i<other._latt_v.size(); i++) {
			_latt_v.push_back(new UnitVisible(*other._latt_v[i]));
		};
		for (auto &pr: other._latt_h) {
			for (auto i=0; i<pr.second.size(); i++) {
				_latt_h[pr.first].push_back(new UnitHidden(*pr.second[i]));
			};
		};
		_latt_v_idxs = other._latt_v_idxs;
		_latt_h_idxs = other._latt_h_idxs;
		for (auto i=0; i<other._conns_hh.size(); i++) {
			_conns_hh.push_back(new ConnHH(*other._conns_hh[i]));
		};
		for (auto i=0; i<other._conns_vv.size(); i++) {
			_conns_vv.push_back(new ConnVV(*other._conns_vv[i]));
		};
		for (auto i=0; i<other._conns_vvv.size(); i++) {
			_conns_vvv.push_back(new ConnVVV(*other._conns_vvv[i]));
		};
		for (auto i=0; i<other._conns_vv.size(); i++) {
			_conns_vh.push_back(new ConnVH(*other._conns_vh[i]));
		};
	};
	void Lattice::_move(Lattice& other) {
		_no_dims = other._no_dims;
		_box_length = other._box_length;
		_latt_v = other._latt_v;
		_latt_h = other._latt_h;
		_conns_hh = other._conns_hh;
		_conns_vv = other._conns_vv;
		_conns_vvv = other._conns_vvv;
		_conns_vh = other._conns_vh;
		_latt_v_idxs = other._latt_v_idxs;
		_latt_h_idxs = other._latt_h_idxs;

		// Reset other
		other._no_dims = 0;
		other._box_length = 0;
		other._latt_v.clear();
		other._latt_h.clear();
		other._conns_hh.clear();
		other._conns_vv.clear();
		other._conns_vvv.clear();
		other._conns_vh.clear();
		other._latt_v_idxs.clear();
		other._latt_h_idxs.clear();
	};

	/********************
	Check setup
	********************/

	void Lattice::print() const {
		// Go through all sites
		for (auto const &s: _latt_v) {
			s->print();
		};
		for (auto const &pr: _latt_h) {
			std::cout << "Layer: " << pr.first << std::endl;
			for (auto const &s: pr.second) {
				s->print();
			};
		};
	};

	/********************
	Helpers to setup all sites - Add possible species to all sites
	********************/

	void Lattice::all_units_v_add_possible_species(Sptr species) {
		if (!species) {
			std::cerr << ">>> Error: Lattice::all_units_v_add_possible_species <<< nullptr is not allowed!" << std::endl;
			exit(EXIT_FAILURE);
		};

		for (auto &s: _latt_v) {
			s->add_possible_species(species);
		};
	};

	void Lattice::all_units_h_add_possible_species(Sptr species) {
		if (!species) {
			std::cerr << ">>> Error: Lattice::all_units_h_add_possible_species <<< nullptr is not allowed!" << std::endl;
			exit(EXIT_FAILURE);
		};

		for (auto &pr: _latt_h) {
			for (auto &s: pr.second) {
				s->add_possible_species(species);
			};
		};	
	};

	/********************
	Helpers to setup all sites - Biases
	********************/

	void Lattice::all_units_v_set_bias_dict(std::shared_ptr<BiasDict> bias_dict) {
		for (auto &s: _latt_v) {
			s->set_bias_dict(bias_dict);
		};
	};

	void Lattice::all_units_h_set_bias_dict(std::shared_ptr<BiasDict> bias_dict) {
		for (auto &pr: _latt_h) {
			for (auto &s: pr.second) {
				s->set_bias_dict(bias_dict);
			};
		};
	};

	/********************
	Helpers to setup all sites - Visible-Visible ixns
	********************/

	void Lattice::all_conns_vv_init() {
		all_conns_vv_init(nullptr);
	};
	void Lattice::all_conns_vv_init(std::shared_ptr<O2IxnDict> ixn_dict) {
		// Neighbor
		UnitVisible *nbr = nullptr;

		// Lattice sites
		for (auto const &s: _latt_v) {

			// Only connect to "plus one" (not minus one) => no duplicates!

			// Dim 1,2,3
			if (s->x()+1 <= _box_length) {
				if (_no_dims == 1) {
					nbr = _look_up_unit_v(s->x()+1);
				} else if (_no_dims == 2) {
					nbr = _look_up_unit_v(s->x()+1,s->y());
				} else if (_no_dims == 3) {
					nbr = _look_up_unit_v(s->x()+1,s->y(),s->z());
				};
				add_conn_vv(s,nbr,ixn_dict);
			};

			// Dim 2,3
			if (_no_dims == 2 || _no_dims == 3) {
				if (s->y()+1 <= _box_length) {
					if (_no_dims == 2) {
						nbr = _look_up_unit_v(s->x(),s->y()+1);
					} else if (_no_dims == 3) {
						nbr = _look_up_unit_v(s->x(),s->y()+1,s->z());
					};
					add_conn_vv(s,nbr,ixn_dict);
				};
			};

			// Dim 3
			if (_no_dims == 3) {
				if (s->z()+1 <= _box_length) {
					if (_no_dims == 3) {
						nbr = _look_up_unit_v(s->x(),s->y(),s->z()+1);
					};
					add_conn_vv(s,nbr,ixn_dict);
				};
			};
		};
	};

	void Lattice::all_conns_vvv_init() {
		all_conns_vvv_init(nullptr);
	};

	void Lattice::all_conns_vvv_init(std::shared_ptr<O3IxnDict> ixn_dict) {

		// Neighbors
		UnitVisible *nbr1 = nullptr, *nbr2 = nullptr;

		// Go through all lattice sites
		for (auto const &s: _latt_v) {

			// Only connect to "plus one/two" (not minus one/two) => no duplicates!

			// Dim 1,2,3
			// Right 1, Right 1
			if (s->x()+2 <= _box_length) {
				if (_no_dims == 1) {
					nbr1 = _look_up_unit_v(s->x()+1);
					nbr2 = _look_up_unit_v(s->x()+2);
				} else if (_no_dims == 2) {
					nbr1 = _look_up_unit_v(s->x()+1,s->y());
					nbr2 = _look_up_unit_v(s->x()+2,s->y());
				} else if (_no_dims == 3) {
					nbr1 = _look_up_unit_v(s->x()+1,s->y(),s->z());
					nbr2 = _look_up_unit_v(s->x()+2,s->y(),s->z());
				};
				add_conn_vvv(s,nbr1,nbr2,ixn_dict);
			};

			// Dim 2,3
			if (_no_dims == 2 || _no_dims == 3) {
				// Up 1, up 1
				if (s->y()+2 <= _box_length) {
					if (_no_dims == 2) {
						nbr1 = _look_up_unit_v(s->x(),s->y()+1);
						nbr2 = _look_up_unit_v(s->x(),s->y()+2);
					} else if (_no_dims == 3) {
						nbr1 = _look_up_unit_v(s->x(),s->y()+1,s->z());
						nbr2 = _look_up_unit_v(s->x(),s->y()+2,s->z());
					};
					add_conn_vvv(s,nbr1,nbr2,ixn_dict);
				};
				// Up 1, right 1
				if (s->x()+1 <= _box_length && s->y()+1 <= _box_length) {
					if (_no_dims == 2) {
						nbr1 = _look_up_unit_v(s->x(),s->y()+1);
						nbr2 = _look_up_unit_v(s->x()+1,s->y()+1);
					} else if (_no_dims == 3) {
						nbr1 = _look_up_unit_v(s->x(),s->y()+1,s->z());
						nbr2 = _look_up_unit_v(s->x()+1,s->y()+1,s->z());
					};
					add_conn_vvv(s,nbr1,nbr2,ixn_dict);
				};
				// Right 1, up 1
				if (s->x()+1 <= _box_length && s->y()+1 <= _box_length) {
					if (_no_dims == 2) {
						nbr1 = _look_up_unit_v(s->x()+1,s->y());
						nbr2 = _look_up_unit_v(s->x()+1,s->y()+1);
					} else if (_no_dims == 3) {
						nbr1 = _look_up_unit_v(s->x()+1,s->y(),s->z());
						nbr2 = _look_up_unit_v(s->x()+1,s->y()+1,s->z());
					};
					add_conn_vvv(s,nbr1,nbr2,ixn_dict);
				};
			};

			// Dim 3
			if (_no_dims == 3) {
				// Forward 1, forward 1
				if (s->z()+2 <= _box_length) {
					if (_no_dims == 3) {
						nbr1 = _look_up_unit_v(s->x(),s->y(),s->z()+1);
						nbr2 = _look_up_unit_v(s->x(),s->y(),s->z()+2);
					};
					add_conn_vvv(s,nbr1,nbr2,ixn_dict);
				};
				// Forward 1, right 1
				if (s->x()+1 <= _box_length && s->z()+1 <= _box_length) {
					if (_no_dims == 3) {
						nbr1 = _look_up_unit_v(s->x(),s->y(),s->z()+1);
						nbr2 = _look_up_unit_v(s->x()+1,s->y(),s->z()+1);
					};
					add_conn_vvv(s,nbr1,nbr2,ixn_dict);
				};
				// Right 1, forward 1
				if (s->x()+1 <= _box_length && s->z()+1 <= _box_length) {
					if (_no_dims == 3) {
						nbr1 = _look_up_unit_v(s->x()+1,s->y(),s->z());
						nbr2 = _look_up_unit_v(s->x()+1,s->y(),s->z()+1);
					};
					add_conn_vvv(s,nbr1,nbr2,ixn_dict);
				};
				// Forward 1, up 1
				if (s->y()+1 <= _box_length && s->z()+1 <= _box_length) {
					if (_no_dims == 3) {
						nbr1 = _look_up_unit_v(s->x(),s->y(),s->z()+1);
						nbr2 = _look_up_unit_v(s->x(),s->y()+1,s->z()+1);
					};
					add_conn_vvv(s,nbr1,nbr2,ixn_dict);
				};
				// Up 1, forward 1
				if (s->y()+1 <= _box_length && s->z()+1 <= _box_length) {
					if (_no_dims == 3) {
						nbr1 = _look_up_unit_v(s->x(),s->y()+1,s->z());
						nbr2 = _look_up_unit_v(s->x(),s->y()+1,s->z()+1);
					};
					add_conn_vvv(s,nbr1,nbr2,ixn_dict);
				};
			};
		};
	};

	/********************
	Helpers to setup all sites - Set ixn dicts of connections
	********************/

	void Lattice::all_conns_vv_set_ixn_dict(std::shared_ptr<O2IxnDict> ixn_dict) {
		for (auto &conn: _conns_vv) {
			conn->set_ixn_dict(ixn_dict);
		};
	};
	void Lattice::all_conns_vvv_set_ixn_dict(std::shared_ptr<O3IxnDict> ixn_dict) {
		for (auto &conn: _conns_vvv) {
			conn->set_ixn_dict(ixn_dict);
		};
	};
	void Lattice::all_conns_vh_set_ixn_dict(std::shared_ptr<O2IxnDict> ixn_dict) {
		for (auto &conn: _conns_vh) {
			conn->set_ixn_dict(ixn_dict);
		};
	};

	/********************
	Helpers to setup all sites - Link units to moments
	********************/

	void Lattice::all_units_v_add_to_moment_h(std::shared_ptr<Moment> moment) {
		for (auto &s: _latt_v) {
			moment->add_unit_to_monitor_h(s);
		};
	};
	void Lattice::all_units_h_add_to_moment_b(std::shared_ptr<Moment> moment) {
		for (auto &pr: _latt_h) {
			for (auto &s: pr.second) {
				moment->add_unit_to_monitor_b(s);
			};
		};
	};
	void Lattice::all_conns_vv_add_to_moment_j(std::shared_ptr<Moment> moment) {
		for (auto &c: _conns_vv) {
			moment->add_conn_to_monitor_j(c);
		};
	};
	void Lattice::all_conns_vvv_add_to_moment_k(std::shared_ptr<Moment> moment) {
		for (auto &c: _conns_vvv) {
			moment->add_conn_to_monitor_k(c);
		};
	};
	void Lattice::all_conns_vh_add_to_moment_w(std::shared_ptr<Moment> moment) {
		for (auto &c: _conns_vh) {
			moment->add_conn_to_monitor_w(c);
		};
	};

	/********************
	Add visible-visible connections
	********************/

	ConnVV* Lattice::add_conn_vv(UnitVisible *uv1, UnitVisible *uv2) {
		_conns_vv.push_back(new ConnVV(uv1,uv2));
		uv1->add_conn(_conns_vv.back(),0);
		uv2->add_conn(_conns_vv.back(),1);
		return _conns_vv.back();
	};
	ConnVV* Lattice::add_conn_vv(UnitVisible *uv1, UnitVisible *uv2, std::shared_ptr<O2IxnDict> ixn_dict) {
		add_conn_vv(uv1,uv2);
		if (ixn_dict) {
			_conns_vv.back()->set_ixn_dict(ixn_dict);
		};
		return _conns_vv.back();
	};
	ConnVVV* Lattice::add_conn_vvv(UnitVisible *uv1, UnitVisible *uv2, UnitVisible *uv3) {
		_conns_vvv.push_back(new ConnVVV(uv1,uv2,uv3));
		uv1->add_conn(_conns_vvv.back(),0);
		uv2->add_conn(_conns_vvv.back(),1);
		uv3->add_conn(_conns_vvv.back(),2);
		return _conns_vvv.back();
	};
	ConnVVV* Lattice::add_conn_vvv(UnitVisible *uv1, UnitVisible *uv2, UnitVisible *uv3, std::shared_ptr<O3IxnDict> ixn_dict) {
		add_conn_vvv(uv1,uv2,uv3);
		if (ixn_dict) {
			_conns_vvv.back()->set_ixn_dict(ixn_dict);
		};
		return _conns_vvv.back();
	};

	/********************
	Add hidden units
	********************/

	UnitHidden* Lattice::add_hidden_unit(int layer) {
		return add_hidden_unit(layer,{});
	};
	UnitHidden* Lattice::add_hidden_unit(int layer, std::vector<Sptr> species_possible) {
		auto it = _latt_h.find(layer);
		int idx;
		if (it != _latt_h.end()) {
			idx = it->second.size();
		} else {
			idx = 0;
		};
		_latt_h[layer].push_back(new UnitHidden(layer,idx,species_possible));
		_latt_h_idxs[layer].push_back(_latt_h[layer].size()-1);
		return _latt_h[layer].back();
	};

	/********************
	Add visible-hidden connections
	********************/

	ConnVH* Lattice::add_conn_vh(UnitVisible *uv, UnitHidden *uh) { 
		_conns_vh.push_back(new ConnVH(uv,uh));
		uv->add_conn(_conns_vh.back());
		uh->add_conn(_conns_vh.back());
		return _conns_vh.back();
	};
	ConnVH* Lattice::add_conn_vh(UnitVisible *uv, UnitHidden *uh, std::shared_ptr<O2IxnDict> ixn_dict) {
		add_conn_vh(uv,uh);
		if (ixn_dict) {
			_conns_vh.back()->set_ixn_dict(ixn_dict);
		};
		return _conns_vh.back();
	};

	/********************
	Add hidden-hidden connections
	********************/

	ConnHH* Lattice::add_conn_hh(UnitHidden *uh1, int layer_1, UnitHidden *uh2, int layer_2) {
		_conns_hh.push_back(new ConnHH(uh1,uh2));
		uh1->add_conn(_conns_hh.back(),0,layer_2);
		uh2->add_conn(_conns_hh.back(),1,layer_1);
		return _conns_hh.back();
	};
	ConnHH* Lattice::add_conn_hh(UnitHidden *uh1, int layer_1, UnitHidden *uh2, int layer_2, std::shared_ptr<O2IxnDict> ixn_dict) {
		add_conn_hh(uh1,layer_1,uh2,layer_2);
		if (ixn_dict) {
			_conns_hh.back()->set_ixn_dict(ixn_dict);
		};
		return _conns_hh.back();
	};

	/********************
	Get unit
	********************/

	const std::vector<UnitVisible*>& Lattice::get_all_units_v() const {
		return _latt_v;
	};
	const std::vector<UnitHidden*>& Lattice::get_all_units_h(int layer) const {
		auto it = _latt_h.find(layer);
		if (it != _latt_h.end()) {
			return it->second;
		};

		// Never get here
		std::cerr << ">>> Error: Lattice::get_all_units_h <<< no units in layer: " << layer << std::endl;
		exit(EXIT_FAILURE);
	};
	const std::map<int,std::vector<UnitHidden*>>& Lattice::get_all_units_h() const {
		return _latt_h;
	};

	UnitVisible* Lattice::get_unit_v(int x) const {
		_check_dim(1);
		return _look_up_unit_v(x);
	};
	UnitVisible* Lattice::get_unit_v(int x, int y) const {
		_check_dim(2);
		return _look_up_unit_v(x,y);
	};
	UnitVisible* Lattice::get_unit_v(int x, int y, int z) const {
		_check_dim(3);
		return _look_up_unit_v(x,y,z);
	};

	UnitHidden* Lattice::get_unit_h(int layer, int idx) const {
		return _look_up_unit_h(layer,idx);
	};

	/********************
	Get connection
	********************/

	ConnVV* Lattice::get_conn_vv(int x1, int x2) const {
		_check_dim(1);

		UnitVisible *uv1 = _look_up_unit_v(x1);
		UnitVisible *uv2 = _look_up_unit_v(x2);

		std::vector<std::pair<ConnVV*,int>> conns = uv1->get_conns_vv();
		for (auto &conn: conns) {
			if (conn.first->check_connects_units(uv1,uv2)) {
				return conn.first;
			};
		};

		// Never get here
		std::cerr << ">>> Error: Lattice::get_conn_vv <<< Connection: " << x1 << " to " << x2 << " not found!" << std::endl;
		exit(EXIT_FAILURE);
	};
	ConnVV* Lattice::get_conn_vv(int x1, int y1, int x2, int y2) const {
		_check_dim(2);

		UnitVisible *uv1 = _look_up_unit_v(x1,y1);
		UnitVisible *uv2 = _look_up_unit_v(x2,y2);

		std::vector<std::pair<ConnVV*,int>> conns = uv1->get_conns_vv();
		for (auto &conn: conns) {
			if (conn.first->check_connects_units(uv1,uv2)) {
				return conn.first;
			};
		};

		// Never get here
		std::cerr << ">>> Error: Lattice::get_conn_vv <<< Connection: " << x1 << " " << y1 << " to " << x2 << " " << y2 << " not found!" << std::endl;
		exit(EXIT_FAILURE);
	};
	ConnVV* Lattice::get_conn_vv(int x1, int y1, int z1, int x2, int y2, int z2) const {
		_check_dim(3);

		UnitVisible *uv1 = _look_up_unit_v(x1,y1,z1);
		UnitVisible *uv2 = _look_up_unit_v(x2,y2,z2);

		std::vector<std::pair<ConnVV*,int>> conns = uv1->get_conns_vv();
		for (auto &conn: conns) {
			if (conn.first->check_connects_units(uv1,uv2)) {
				return conn.first;
			};
		};

		// Never get here
		std::cerr << ">>> Error: Lattice::get_conn_vv <<< Connection: " << x1 << " " << y1 << " " << z1 << " to " << x2 << " " << y2 << " " << z2 << " not found!" << std::endl;
		exit(EXIT_FAILURE);
	};

	ConnVVV* Lattice::get_conn_vvv(int x1, int x2, int x3) const {
		_check_dim(1);

		UnitVisible *uv1 = _look_up_unit_v(x1);
		UnitVisible *uv2 = _look_up_unit_v(x2);
		UnitVisible *uv3 = _look_up_unit_v(x3);

		std::vector<std::pair<ConnVVV*,int>> conns = uv1->get_conns_vvv();
		for (auto &conn: conns) {
			if (conn.first->check_connects_units(uv1,uv2,uv3)) {
				return conn.first;
			};
		};

		// Never get here
		std::cerr << ">>> Error: Lattice::get_conn_vvv <<< Connection: " << x1 << " to " << x2 << " to " << x3 << " not found!" << std::endl;
		exit(EXIT_FAILURE);
	};
	ConnVVV* Lattice::get_conn_vvv(int x1, int y1, int x2, int y2, int x3, int y3) const {
		_check_dim(2);

		UnitVisible *uv1 = _look_up_unit_v(x1,y1);
		UnitVisible *uv2 = _look_up_unit_v(x2,y2);
		UnitVisible *uv3 = _look_up_unit_v(x3,y3);

		std::vector<std::pair<ConnVVV*,int>> conns = uv1->get_conns_vvv();
		for (auto &conn: conns) {
			if (conn.first->check_connects_units(uv1,uv2,uv3)) {
				return conn.first;
			};
		};

		// Never get here
		std::cerr << ">>> Error: Lattice::get_conn_vvv <<< Connection: " << x1 << " " << y1 << " to " << x2 << " " << y2 << " to " << x3 << " " << y3 << " not found!" << std::endl;
		exit(EXIT_FAILURE);	
	};
	ConnVVV* Lattice::get_conn_vvv(int x1, int y1, int z1, int x2, int y2, int z2, int x3, int y3, int z3) const {
		_check_dim(3);

		UnitVisible *uv1 = _look_up_unit_v(x1,y1,z1);
		UnitVisible *uv2 = _look_up_unit_v(x2,y2,z2);
		UnitVisible *uv3 = _look_up_unit_v(x3,y3,z3);

		std::vector<std::pair<ConnVVV*,int>> conns = uv1->get_conns_vvv();
		for (auto &conn: conns) {
			if (conn.first->check_connects_units(uv1,uv2,uv3)) {
				return conn.first;
			};
		};

		// Never get here
		std::cerr << ">>> Error: Lattice::get_conn_vvv <<< Connection: " << x1 << " " << y1 << " " << z1 << " to " << x2 << " " << y2 << " " << z2 << " to " << x3 << " " << y3 << " " << z3 << " not found!" << std::endl;
		exit(EXIT_FAILURE);
	};

	ConnVH* Lattice::get_conn_vh(int x, int layer, int idx) const {
		_check_dim(1);

		UnitVisible *uv = _look_up_unit_v(x);
		UnitHidden *uh = _look_up_unit_h(layer,idx);

		std::vector<ConnVH*> conns = uv->get_conns_vh();
		for (auto &conn: conns) {
			if (conn->check_connects_units(uv,uh)) {
				return conn;
			};
		};

		// Never get here
		std::cerr << ">>> Error: Lattice::get_conn_vh <<< Connection: " << x << " to " << layer << " " << idx << " not found!" << std::endl;
		exit(EXIT_FAILURE);
	};
	ConnVH* Lattice::get_conn_vh(int x, int y, int layer, int idx) const {
		_check_dim(2);

		UnitVisible *uv = _look_up_unit_v(x,y);
		UnitHidden *uh = _look_up_unit_h(layer,idx);

		std::vector<ConnVH*> conns = uv->get_conns_vh();
		for (auto &conn: conns) {
			if (conn->check_connects_units(uv,uh)) {
				return conn;
			};
		};

		// Never get here
		std::cerr << ">>> Error: Lattice::get_conn_vh <<< Connection: " << x << " " << y << " to " << layer << " " << idx << " not found!" << std::endl;
		exit(EXIT_FAILURE);
	};
	ConnVH* Lattice::get_conn_vh(int x, int y, int z, int layer, int idx) const {
		_check_dim(3);

		UnitVisible *uv = _look_up_unit_v(x,y,z);
		UnitHidden *uh = _look_up_unit_h(layer,idx);

		std::vector<ConnVH*> conns = uv->get_conns_vh();
		for (auto &conn: conns) {
			if (conn->check_connects_units(uv,uh)) {
				return conn;
			};
		};

		// Never get here
		std::cerr << ">>> Error: Lattice::get_conn_vh <<< Connection: " << x << " " << y << " " << z << " to " << layer << " " << idx << " not found!" << std::endl;
		exit(EXIT_FAILURE);
	};

	ConnHH* Lattice::get_conn_hh(int layer1, int idx1, int layer2, int idx2) const {
		UnitHidden *uh1 = _look_up_unit_h(layer1,idx1);
		UnitHidden *uh2 = _look_up_unit_h(layer2,idx2);

		std::map<int,std::vector<std::pair<ConnHH*,int>>> conns = uh1->get_conns_hh();
		auto it = conns.find(layer1);
		if (it != conns.end()) {
			for (auto &conn: it->second) {
				if (conn.first->check_connects_units(uh1,uh2)) {
					return conn.first;
				};
			};
		};

		// Never get here
		std::cerr << ">>> Error: Lattice::get_conn_hh <<< Connection: " << layer1 << " " << idx1 << " to " << layer2 << " " << idx2 << " not found!" << std::endl;
		exit(EXIT_FAILURE);
	};

	/********************
	Getters
	********************/

	int Lattice::get_no_dims() const {
		return _no_dims;
	};
	int Lattice::get_no_units_v() { 
		return _latt_v.size(); 
	};
	int Lattice::get_no_units_h() { 
		int count=0;
		for (auto const &pr: _latt_h) {
			count += pr.second.size();
		};
		return count; 
	};

	/********************
	Apply funcs to all units
	********************/

	// Clear the lattice
	void Lattice::all_units_v_set_empty() {
		for (auto &s: _latt_v) {
			s->set_b_mode_empty();
			s->set_p_mode_empty();
		};
	};
	void Lattice::all_units_h_set_empty() {
		for (auto &pr: _latt_h) {
			for (auto &s: pr.second) {
				s->set_b_mode_empty();
				s->set_p_mode_empty();
			};
		};
	};
	void Lattice::all_units_set_empty() {
		all_units_v_set_empty();
		all_units_h_set_empty();
	};

	// Random
	void Lattice::all_units_v_random() {
		all_units_set_empty();

		for (auto &ptr: _latt_v) {
			ptr->set_b_mode_random();
		};
	};

	// Binary/probabilistic
	void Lattice::all_units_convert_to_b_mode() {
		for (auto &s: _latt_v) {
			s->convert_p_to_b_mode();
		};
		for (auto &pr: _latt_h) {
			for (auto &s: pr.second) {
				s->convert_p_to_b_mode();
			};
		};
	};
	void Lattice::all_units_convert_to_p_mode() {
		for (auto &s: _latt_v) {
			s->convert_b_to_p_mode();
		};
		for (auto &pr: _latt_h) {
			for (auto &s: pr.second) {
				s->convert_b_to_p_mode();
			};
		};
	};

	/********************
	Write/read Latt to a file
	********************/

	void Lattice::write_to_file(std::string fname, bool binary) 
	{
		std::ofstream f;
		f.open (fname);
		for (auto const &l: _latt_v) {
			if (binary) {
				if (!l->check_b_mode_is_empty()) {
					if (_no_dims == 1) {
						f << l->x();
					} else if (_no_dims == 2) {
						f << l->x() << " " << l->y();
					} else if (_no_dims == 3) {
						f << l->x() << " " << l->y() << " " << l->z();
					};
					f << " " << l->get_b_mode_species()->get_name() << "\n";
				};
			} else {
				if (!l->check_p_mode_is_empty()) {
					if (_no_dims == 1) {
						f << l->x();
					} else if (_no_dims == 2) {
						f << l->x() << " " << l->y();
					} else if (_no_dims == 3) {
						f << l->x() << " " << l->y() << " " << l->z();
					};
					for (auto const &pr: l->get_p_mode_probs()) {
						f << " " << pr.first->get_name() << " " << pr.second;
					};
					f << "\n";
				};
			};
		};
		f.close();
	};
	void Lattice::read_from_file(std::string fname, bool binary)
	{
		// Clear the current lattice
		all_units_set_empty();

		std::ifstream f;
		f.open(fname);
		std::string x="",y="",z="";
		std::string sp="";
		std::string line;
		std::istringstream iss;
		UnitVisible* s;
		std::string prob="";
		double prob_val;
		if (f.is_open()) { // make sure we found it
			while (getline(f,line)) {
				if (line == "") { continue; };
				iss = std::istringstream(line);
			    iss >> x;
			    if (_no_dims == 2 || _no_dims == 3) {
				    iss >> y;
			    };
			    if (_no_dims == 3) {
				    iss >> z;
			    };
			    iss >> sp;
			    if (!binary) {
			    	// Read the prob
			    	iss >> prob;
			    };
		    	// Add to Latt
		    	if (_no_dims == 1) {
		    		s = _look_up_unit_v(atoi(x.c_str()));
			    } else if (_no_dims == 2) {
		    		s = _look_up_unit_v(atoi(x.c_str()),atoi(y.c_str()));
			    } else if (_no_dims == 3) {
			    	s = _look_up_unit_v(atoi(x.c_str()),atoi(y.c_str()),atoi(z.c_str()));
			    };
			    if (binary) {
		    		s->set_b_mode_species(sp);
		    	} else {
		    		prob_val = atof(prob.c_str());
		    		s->set_p_mode_prob(sp,prob_val);
		    		s->set_p_mode_prob("",1.0-prob_val);
		    	};
	    		// Reset
		    	sp=""; x=""; y=""; z=""; prob="";
			};
		};
		f.close();

		// std::cout << "Read: " << fname << std::endl;
	};

	/********************
	Sample
	********************/

	void Lattice::sample_v(bool binary) {
		
		// Go through hidden layers top to bottom (reverse)
		auto rit_begin = _latt_h.rbegin();
		rit_begin++; // Skip the first
		for (auto rit=rit_begin; rit!=_latt_h.rend(); ++rit)
		{
			// Shuffle
			std::random_shuffle ( _latt_h_idxs[rit->first].begin(), _latt_h_idxs[rit->first].end() );

			for (auto const &idx: _latt_h_idxs[rit->first]) 
			{
				// Sample, given input from a layer higher
				rit->second[idx]->sample(rit->first+1, binary);
			};
		};

		// Finally, sample the visible layer

		// Shuffle
		std::random_shuffle ( _latt_v_idxs.begin(), _latt_v_idxs.end() );		

		// Sample
		for (auto const &idx: _latt_v_idxs) 
		{
			_latt_v[idx]->sample(binary);
		};
	};
	void Lattice::sample_h(bool binary) {
		for (auto it=_latt_h.begin(); it!=_latt_h.end(); it++) 
		{
			// Shuffle
			std::random_shuffle ( _latt_h_idxs[it->first].begin(), _latt_h_idxs[it->first].end() );

			for (auto const &idx: _latt_h_idxs[it->first]) 
			{
				// Sample, given input from a layer lower
				it->second[idx]->sample(it->first-1,binary);
			};
		};
	};

	/********************
	Get counts
	********************/

	// 1 particle
	double Lattice::get_count(Sptr &sp, bool binary) const {
		double count = 0.0;
		if (binary) {
			for (auto const &s: _latt_v) {
				if (sp == s->get_b_mode_species()) {
					count += 1.0;
				};
			};
		} else {
			for (auto const &s: _latt_v) {
				count += s->get_p_mode_prob(sp);
			};
		};

		return count;
	};

	// 2 particle
	double Lattice::get_count(Sptr &sp1, Sptr &sp2, bool binary, bool reversibly) const {
		// Only dim=1 for now
		if (_no_dims != 1) {
			std::cerr << ">>> Error: Lattice::get_count <<< only supported for dim 1." << std::endl;
			exit(EXIT_FAILURE);
		};

		UnitVisible *nbr = nullptr;
		double count = 0.0;
		for (auto &s: _latt_v) {
			// Only connect to "plus one" (not minus one) => no duplicates!

			if (s->x()+1 <= _box_length) {
				// Get neighbor
				nbr = _look_up_unit_v(s->x()+1);

				// Count
				if (binary) {
					if ((sp1 == s->get_b_mode_species()) && (sp2 == nbr->get_b_mode_species())) {
						count += 1.0;
					};
					if (reversibly && sp1 != sp2) {
						if ((sp1 == nbr->get_b_mode_species()) && (sp2 == s->get_b_mode_species())) {
							count += 1.0;
						};
					};
				} else {
					count += s->get_p_mode_prob(sp1) * nbr->get_p_mode_prob(sp2); 
					if (reversibly && sp1 != sp2) {
						count += nbr->get_p_mode_prob(sp1) * s->get_p_mode_prob(sp2); 
					};
				};
			};
		};

		return count;
	};

	// 3 particle
	double Lattice::get_count(Sptr &sp1, Sptr &sp2, Sptr &sp3, bool binary, bool reversibly) const {
		// Only dim=1 for now
		if (_no_dims != 1) {
			std::cerr << ">>> Error: Lattice::get_count <<< only supported for dim 1." << std::endl;
			exit(EXIT_FAILURE);
		};

		UnitVisible *nbr1 = nullptr, *nbr2 = nullptr;
		double count = 0.0;
		for (auto &s: _latt_v) {
			// Only connect to "plus one" (not minus one) => no duplicates!

			if (s->x()+2 <= _box_length) {
				// Get neighbors
				nbr1 = _look_up_unit_v(s->x()+1);
				nbr2 = _look_up_unit_v(s->x()+2);

				// Count
				if (binary) {
					if ((sp1 == s->get_b_mode_species()) && (sp2 == nbr1->get_b_mode_species()) && (sp3 == nbr2->get_b_mode_species())) {
						count += 1.0;
					};
					if (reversibly && sp1 != sp3) {
						if ((sp1 == nbr2->get_b_mode_species()) && (sp2 == nbr1->get_b_mode_species()) && (sp3 == s->get_b_mode_species())) {
							count += 1.0;
						};
					};
				} else {
					count += s->get_p_mode_prob(sp1) * nbr1->get_p_mode_prob(sp2) * nbr2->get_p_mode_prob(sp3); 
					if (reversibly && sp1 != sp3) {
						count += nbr2->get_p_mode_prob(sp1) * nbr1->get_p_mode_prob(sp2) * s->get_p_mode_prob(sp3); 
					};
				};
			};
		};

		return count;
	};

	// 4 particle
	double Lattice::get_count(Sptr &sp1, Sptr &sp2, Sptr &sp3, Sptr &sp4, bool binary, bool reversibly) const {
		// Only dim=1 for now
		if (_no_dims != 1) {
			std::cerr << ">>> Error: Lattice::get_count <<< only supported for dim 1." << std::endl;
			exit(EXIT_FAILURE);
		};

		UnitVisible *nbr1 = nullptr, *nbr2 = nullptr, *nbr3 = nullptr;
		double count = 0.0;
		for (auto &s: _latt_v) {
			// Only connect to "plus one" (not minus one) => no duplicates!

			if (s->x()+3 <= _box_length) {
				// Get neighbors
				nbr1 = _look_up_unit_v(s->x()+1);
				nbr2 = _look_up_unit_v(s->x()+2);
				nbr3 = _look_up_unit_v(s->x()+3);

				// Count
				if (binary) {
					if ((sp1 == s->get_b_mode_species()) && (sp2 == nbr1->get_b_mode_species()) && (sp3 == nbr2->get_b_mode_species()) && (sp4 == nbr3->get_b_mode_species())) {
						count += 1.0;
					};
					if (reversibly && !(sp1 == sp4 && sp2 == sp3)) {
						if ((sp1 == nbr3->get_b_mode_species()) && (sp2 == nbr2->get_b_mode_species()) && (sp3 == nbr1->get_b_mode_species()) && (sp4 == s->get_b_mode_species())) {
							count += 1.0;
						};
					};
				} else {
					count += s->get_p_mode_prob(sp1) * nbr1->get_p_mode_prob(sp2) * nbr2->get_p_mode_prob(sp3) * nbr3->get_p_mode_prob(sp4); 
					if (reversibly && !(sp1 == sp4 && sp2 == sp3)) {
						count += nbr3->get_p_mode_prob(sp1) * nbr2->get_p_mode_prob(sp2) * nbr1->get_p_mode_prob(sp3) * s->get_p_mode_prob(sp4); 
					};
				};
			};
		};

		return count;
	};

	/****************************************
	Lattice - PRIVATE
	****************************************/

	/********************
	Lookup a site iterator from x,y,z
	********************/

	UnitVisible* Lattice::_look_up_unit_v(int x) const {
		if (x > _box_length || x < 1) {
			return nullptr;
		};

		// Figure out index in list
		int n = x-1;

		return _latt_v[n];
	};
	UnitVisible* Lattice::_look_up_unit_v(int x, int y) const {
		if (x > _box_length || x < 1 || y > _box_length || y < 1) {
			return nullptr;
		};

		// Figure out index in list
		int n = (x-1)*_box_length + y-1;

		return _latt_v[n];
	};
	UnitVisible* Lattice::_look_up_unit_v(int x, int y, int z) const {
		if (x > _box_length || x < 1 || y > _box_length || y < 1 || z > _box_length || z < 1) {
			return nullptr;
		};

		// Figure out index in list
		int n = (x-1)*_box_length*_box_length + (y-1)*_box_length + z-1;

		return _latt_v[n];
	};

	UnitHidden* Lattice::_look_up_unit_h(int layer, int idx) const {
		auto it = _latt_h.find(layer);
		if (it == _latt_h.end()) {
			return nullptr;
		};
		if (idx >= it->second.size()) {
			return nullptr;
		};

		return it->second[idx];
	};

	/********************
	Check dim
	********************/

	void Lattice::_check_dim(int dim) const {
		if (_no_dims != dim) {
			std::cerr << ">>> Error: Lattice::_check_dim <<< dim is: " << _no_dims << " but requested: " << dim << std::endl;
			exit(EXIT_FAILURE);
		};	
	};
};