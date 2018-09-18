#include "../../include/dynamicboltz_bits/lattice.hpp"

// Other headers
#include "../../include/dynamicboltz_bits/general.hpp"
#include "../../include/dynamicboltz_bits/species.hpp"
#include "../../include/dynamicboltz_bits/unit.hpp"
#include "../../include/dynamicboltz_bits/connections.hpp"
#include "../../include/dynamicboltz_bits/ixn_dicts.hpp"
#include "../../include/dynamicboltz_bits/moment.hpp"

#include <iostream>
#include <fstream>
#include <numeric>
#include "math.h"
#include <ctime>
#include <sstream>
#include <random>

/************************************
* Namespace for dblz
************************************/

namespace dblz {

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
		_dim = dim;
		_box_length = box_length;

		// Make a fully linked list of sites
		if (dim == 1) {
			for (int x=1; x<=box_length; x++) {
				_latt.push_back(UnitVisible(x));					
			};
		} else if (dim == 2) {
			for (int x=1; x<=box_length; x++) {
				for (int y=1; y<=box_length; y++) {
					_latt.push_back(UnitVisible(x,y));					
				};
			};
		} else if (dim == 3) {
			for (int x=1; x<=box_length; x++) {
				for (int y=1; y<=box_length; y++) {
					for (int z=1; z<=box_length; z++) {
						_latt.push_back(UnitVisible(x,y,z));					
					};
				};
			};
		};
	};
	Lattice::Lattice(int dim, int box_length, std::vector<Sptr> possible_species_of_all_units_vis) : Lattice(dim,box_length) {
		for (Latt_it lit=_latt.begin(); lit != _latt.end(); lit++) {
			lit->set_possible_species(possible_species_of_all_units_vis);
		};
	};
	Lattice::Lattice(const Lattice& other) {
		_copy(other);
	};
	Lattice::Lattice(Lattice&& other) {
		_copy(other);
		other._reset();
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
			_copy(other);
			other._reset();
		};
		return *this;
	};
	Lattice::~Lattice() {
		_clean_up();
	};

	void Lattice::_clean_up() {
		// Nothing...
	};
	void Lattice::_copy(const Lattice& other) {
		_dim = other._dim;
		_box_length = other._box_length;
		_latt = other._latt;
		_conns_vv = other._conns_vv;
		_conns_vvv = other._conns_vvv;
	};
	void Lattice::_reset() {
		_dim = 0;
		_box_length = 0;
		_latt.clear();
		_conns_vv.clear();
		_conns_vvv.clear();
	};

	/********************
	Check setup
	********************/

	void Lattice::check_setup() const {
		// Go through all sites
		for (auto const &s: _latt) {
			s.check_setup();
		};
	};

	void Lattice::print_occupancy(bool binary) const {
		for (auto const &s: _latt) {
			if (_dim == 1) {
				std::cout << s.x() << " " << std::flush;
			} else if (_dim == 2) {
				std::cout << s.x() << " " << s.y() << " " << std::flush;
			} else if (_dim == 3) {
				std::cout << s.x() << " " << s.y() << " " << s.z() << " " << std::flush;
			};
			if (binary) {
				if (s.get_b_mode_species()) {
					std::cout << s.get_b_mode_species()->get_name() << std::endl;
				} else {
					std::cout << "EMPTY" << std::endl;
				};
			} else {
				for (auto &pr: s.get_p_mode_probs()) {
					std::cout << "(" << pr.first->get_name() << " " << pr.second << ") ";
				};
				std::cout << std::endl;
			};
		};
	};

	/********************
	Finish setup
	********************/

	// Add possible species to all sites
	void Lattice::add_possible_species_to_all_units_vis(Sptr species) {
		if (!species) {
			std::cerr << ">>> Error: Lattice::add_possible_species_to_all_units_vis <<< nullptr is not allowed!" << std::endl;
			exit(EXIT_FAILURE);
		};

		for (auto &s: _latt) {
			s.add_possible_species(species);
		};
	};

	// Biases
	void Lattice::set_bias_dict_of_all_units_vis(std::shared_ptr<BiasDict> bias_dict) {
		for (auto &s: _latt) {
			s.set_bias_dict(bias_dict);
		};
	};

	// Visible-Visible ixns
	void Lattice::init_conns_NN_all_units_vis() {
		init_conns_NN_all_units_vis(nullptr);
	};
	void Lattice::init_conns_NN_all_units_vis(std::shared_ptr<O2IxnDict> ixn_dict) {

		Latt_it lit = _latt.begin();
		UnitVisible *nbr = nullptr;
		while (lit != _latt.end()) {
			// Only connect to "plus one" (not minus one) => no duplicates!

			// Dim 1,2,3
			if (lit->x()+1 <= _box_length) {
				if (_dim == 1) {
					nbr = &_look_up(lit->x()+1);
				} else if (_dim == 2) {
					nbr = &_look_up(lit->x()+1,lit->y());
				} else if (_dim == 3) {
					nbr = &_look_up(lit->x()+1,lit->y(),lit->z());
				};
				add_conn_vv(&*lit,nbr,ixn_dict);
			};

			// Dim 2,3
			if (_dim == 2 || _dim == 3) {
				if (lit->y()+1 <= _box_length) {
					if (_dim == 2) {
						nbr = &_look_up(lit->x(),lit->y()+1);
					} else if (_dim == 3) {
						nbr = &_look_up(lit->x(),lit->y()+1,lit->z());
					};
					add_conn_vv(&*lit,nbr,ixn_dict);
				};
			};

			// Dim 3
			if (_dim == 3) {
				if (lit->z()+1 <= _box_length) {
					if (_dim == 3) {
						nbr = &_look_up(lit->x(),lit->y(),lit->z()+1);
					};
					add_conn_vv(&*lit,nbr,ixn_dict);
				};
			};

			// Next
			lit++;
		};
	};

	void Lattice::init_conns_triplet_all_units_vis() {
		init_conns_triplet_all_units_vis(nullptr);
	};
	void Lattice::init_conns_triplet_all_units_vis(std::shared_ptr<O3IxnDict> ixn_dict) {

		Latt_it lit = _latt.begin();
		UnitVisible *nbr1 = nullptr, *nbr2 = nullptr;
		while (lit != _latt.end()) {
			// Only connect to "plus one/two" (not minus one/two) => no duplicates!

			// Dim 1,2,3
			// Right 1, Right 1
			if (lit->x()+2 <= _box_length) {
				if (_dim == 1) {
					nbr1 = &_look_up(lit->x()+1);
					nbr2 = &_look_up(lit->x()+2);
				} else if (_dim == 2) {
					nbr1 = &_look_up(lit->x()+1,lit->y());
					nbr2 = &_look_up(lit->x()+2,lit->y());
				} else if (_dim == 3) {
					nbr1 = &_look_up(lit->x()+1,lit->y(),lit->z());
					nbr2 = &_look_up(lit->x()+2,lit->y(),lit->z());
				};
				add_conn_vvv(&*lit,nbr1,nbr2,ixn_dict);
			};

			// Dim 2,3
			if (_dim == 2 || _dim == 3) {
				// Up 1, up 1
				if (lit->y()+2 <= _box_length) {
					if (_dim == 2) {
						nbr1 = &_look_up(lit->x(),lit->y()+1);
						nbr2 = &_look_up(lit->x(),lit->y()+2);
					} else if (_dim == 3) {
						nbr1 = &_look_up(lit->x(),lit->y()+1,lit->z());
						nbr2 = &_look_up(lit->x(),lit->y()+2,lit->z());
					};
					add_conn_vvv(&*lit,nbr1,nbr2,ixn_dict);
				};
				// Up 1, right 1
				if (lit->x()+1 <= _box_length && lit->y()+1 <= _box_length) {
					if (_dim == 2) {
						nbr1 = &_look_up(lit->x(),lit->y()+1);
						nbr2 = &_look_up(lit->x()+1,lit->y()+1);
					} else if (_dim == 3) {
						nbr1 = &_look_up(lit->x(),lit->y()+1,lit->z());
						nbr2 = &_look_up(lit->x()+1,lit->y()+1,lit->z());
					};
					add_conn_vvv(&*lit,nbr1,nbr2,ixn_dict);
				};
				// Right 1, up 1
				if (lit->x()+1 <= _box_length && lit->y()+1 <= _box_length) {
					if (_dim == 2) {
						nbr1 = &_look_up(lit->x()+1,lit->y());
						nbr2 = &_look_up(lit->x()+1,lit->y()+1);
					} else if (_dim == 3) {
						nbr1 = &_look_up(lit->x()+1,lit->y(),lit->z());
						nbr2 = &_look_up(lit->x()+1,lit->y()+1,lit->z());
					};
					add_conn_vvv(&*lit,nbr1,nbr2,ixn_dict);
				};
			};

			// Dim 3
			if (_dim == 3) {
				// Forward 1, forward 1
				if (lit->z()+2 <= _box_length) {
					if (_dim == 3) {
						nbr1 = &_look_up(lit->x(),lit->y(),lit->z()+1);
						nbr2 = &_look_up(lit->x(),lit->y(),lit->z()+2);
					};
					add_conn_vvv(&*lit,nbr1,nbr2,ixn_dict);
				};
				// Forward 1, right 1
				if (lit->x()+1 <= _box_length && lit->z()+1 <= _box_length) {
					if (_dim == 3) {
						nbr1 = &_look_up(lit->x(),lit->y(),lit->z()+1);
						nbr2 = &_look_up(lit->x()+1,lit->y(),lit->z()+1);
					};
					add_conn_vvv(&*lit,nbr1,nbr2,ixn_dict);
				};
				// Right 1, forward 1
				if (lit->x()+1 <= _box_length && lit->z()+1 <= _box_length) {
					if (_dim == 3) {
						nbr1 = &_look_up(lit->x()+1,lit->y(),lit->z());
						nbr2 = &_look_up(lit->x()+1,lit->y(),lit->z()+1);
					};
					add_conn_vvv(&*lit,nbr1,nbr2,ixn_dict);
				};
				// Forward 1, up 1
				if (lit->y()+1 <= _box_length && lit->z()+1 <= _box_length) {
					if (_dim == 3) {
						nbr1 = &_look_up(lit->x(),lit->y(),lit->z()+1);
						nbr2 = &_look_up(lit->x(),lit->y()+1,lit->z()+1);
					};
					add_conn_vvv(&*lit,nbr1,nbr2,ixn_dict);
				};
				// Up 1, forward 1
				if (lit->y()+1 <= _box_length && lit->z()+1 <= _box_length) {
					if (_dim == 3) {
						nbr1 = &_look_up(lit->x(),lit->y()+1,lit->z());
						nbr2 = &_look_up(lit->x(),lit->y()+1,lit->z()+1);
					};
					add_conn_vvv(&*lit,nbr1,nbr2,ixn_dict);
				};
			};

			// Next
			lit++;
		};
	};

	void Lattice::set_ixn_dict_of_all_conns_vv(std::shared_ptr<O2IxnDict> ixn_dict) {
		for (auto &conn: _conns_vv) {
			conn.set_ixn_dict(ixn_dict);
		};
	};
	void Lattice::set_ixn_dict_of_all_conns_vvv(std::shared_ptr<O3IxnDict> ixn_dict) {
		for (auto &conn: _conns_vvv) {
			conn.set_ixn_dict(ixn_dict);
		};
	};

	void Lattice::add_conn_vv(UnitVisible *uv1, UnitVisible *uv2, std::shared_ptr<O2IxnDict> ixn_dict) {
		_conns_vv.push_back(ConnVV(uv1,uv2));
		if (ixn_dict) {
			_conns_vv.back().set_ixn_dict(ixn_dict);
		};
		uv1->add_conn(&_conns_vv.back(),0);
		uv2->add_conn(&_conns_vv.back(),1);
	};
	void Lattice::add_conn_vvv(UnitVisible *uv1, UnitVisible *uv2, UnitVisible *uv3, std::shared_ptr<O3IxnDict> ixn_dict) {
		_conns_vvv.push_back(ConnVVV(uv1,uv2,uv3));
		if (ixn_dict) {
			_conns_vvv.back().set_ixn_dict(ixn_dict);
		};
		uv1->add_conn(&_conns_vvv.back(),0);
		uv2->add_conn(&_conns_vvv.back(),1);
		uv3->add_conn(&_conns_vvv.back(),2);
	};

	void Lattice::add_all_units_vis_to_moment_h(std::shared_ptr<Moment> moment) {
		for (auto &s: _latt) {
			moment->add_unit_to_monitor_h(&s);
		};
	};
	void Lattice::add_all_conns_vv_to_moment_j(std::shared_ptr<Moment> moment) {
		for (auto &c: _conns_vv) {
			moment->add_conn_to_monitor_j(&c);
		};
	};
	void Lattice::add_all_conns_vvv_to_moment_k(std::shared_ptr<Moment> moment) {
		for (auto &c: _conns_vvv) {
			moment->add_conn_to_monitor_k(&c);
		};
	};

	/********************
	Get unit
	********************/

	UnitVisible& Lattice::get_unit_visible(int x) {
		_check_dim(1);
		return _look_up(x);
	};
	UnitVisible& Lattice::get_unit_visible(int x, int y) {
		_check_dim(2);
		return _look_up(x,y);
	};
	UnitVisible& Lattice::get_unit_visible(int x, int y, int z) {
		_check_dim(3);
		return _look_up(x,y,z);
	};

	/********************
	Get connection
	********************/

	ConnVV& Lattice::get_conn_vv(int x1, int x2) {
		_check_dim(1);

		UnitVisible *uv1 = &_look_up(x1);
		UnitVisible *uv2 = &_look_up(x2);

		std::vector<std::pair<ConnVV*,int>> conns = uv1->get_conns_vv();
		for (auto &conn: conns) {
			if (conn.first->check_connects_units(uv1,uv2)) {
				return *(conn.first);
			};
		};

		// Never get here
		std::cerr << ">>> Error: Lattice::get_conn_vv <<< Connection: " << x1 << " to " << x2 << " not found!" << std::endl;
		exit(EXIT_FAILURE);
	};
	ConnVV& Lattice::get_conn_vv(int x1, int y1, int x2, int y2) {
		_check_dim(2);

		UnitVisible *uv1 = &_look_up(x1,y1);
		UnitVisible *uv2 = &_look_up(x2,y2);

		std::vector<std::pair<ConnVV*,int>> conns = uv1->get_conns_vv();
		for (auto &conn: conns) {
			if (conn.first->check_connects_units(uv1,uv2)) {
				return *(conn.first);
			};
		};

		// Never get here
		std::cerr << ">>> Error: Lattice::get_conn_vv <<< Connection: " << x1 << " " << y1 << " to " << x2 << " " << y2 << " not found!" << std::endl;
		exit(EXIT_FAILURE);
	};
	ConnVV& Lattice::get_conn_vv(int x1, int y1, int z1, int x2, int y2, int z2) {
		_check_dim(3);

		UnitVisible *uv1 = &_look_up(x1,y1,z1);
		UnitVisible *uv2 = &_look_up(x2,y2,z2);

		std::vector<std::pair<ConnVV*,int>> conns = uv1->get_conns_vv();
		for (auto &conn: conns) {
			if (conn.first->check_connects_units(uv1,uv2)) {
				return *(conn.first);
			};
		};

		// Never get here
		std::cerr << ">>> Error: Lattice::get_conn_vv <<< Connection: " << x1 << " " << y1 << " " << z1 << " to " << x2 << " " << y2 << " " << z2 << " not found!" << std::endl;
		exit(EXIT_FAILURE);
	};

	ConnVVV& Lattice::get_conn_vvv(int x1, int x2, int x3) {
		_check_dim(1);

		UnitVisible *uv1 = &_look_up(x1);
		UnitVisible *uv2 = &_look_up(x2);
		UnitVisible *uv3 = &_look_up(x3);

		std::vector<std::pair<ConnVVV*,int>> conns = uv1->get_conns_vvv();
		for (auto &conn: conns) {
			if (conn.first->check_connects_units(uv1,uv2,uv3)) {
				return *(conn.first);
			};
		};

		// Never get here
		std::cerr << ">>> Error: Lattice::get_conn_vvv <<< Connection: " << x1 << " to " << x2 << " to " << x3 << " not found!" << std::endl;
		exit(EXIT_FAILURE);
	};
	ConnVVV& Lattice::get_conn_vvv(int x1, int y1, int x2, int y2, int x3, int y3) {
		_check_dim(2);

		UnitVisible *uv1 = &_look_up(x1,y1);
		UnitVisible *uv2 = &_look_up(x2,y2);
		UnitVisible *uv3 = &_look_up(x3,y3);

		std::vector<std::pair<ConnVVV*,int>> conns = uv1->get_conns_vvv();
		for (auto &conn: conns) {
			if (conn.first->check_connects_units(uv1,uv2,uv3)) {
				return *(conn.first);
			};
		};

		// Never get here
		std::cerr << ">>> Error: Lattice::get_conn_vvv <<< Connection: " << x1 << " " << y1 << " to " << x2 << " " << y2 << " to " << x3 << " " << y3 << " not found!" << std::endl;
		exit(EXIT_FAILURE);	
	};
	ConnVVV& Lattice::get_conn_vvv(int x1, int y1, int z1, int x2, int y2, int z2, int x3, int y3, int z3) {
		_check_dim(3);

		UnitVisible *uv1 = &_look_up(x1,y1,z1);
		UnitVisible *uv2 = &_look_up(x2,y2,z2);
		UnitVisible *uv3 = &_look_up(x3,y3,z3);

		std::vector<std::pair<ConnVVV*,int>> conns = uv1->get_conns_vvv();
		for (auto &conn: conns) {
			if (conn.first->check_connects_units(uv1,uv2,uv3)) {
				return *(conn.first);
			};
		};

		// Never get here
		std::cerr << ">>> Error: Lattice::get_conn_vvv <<< Connection: " << x1 << " " << y1 << " " << z1 << " to " << x2 << " " << y2 << " " << z2 << " to " << x3 << " " << y3 << " " << z3 << " not found!" << std::endl;
		exit(EXIT_FAILURE);
	};

	/********************
	Getters
	********************/

	int Lattice::dim() const {
		return _dim;
	};
	int Lattice::no_units_vis() { 
		return _latt.size(); 
	};

	/********************
	Apply funcs to all units
	********************/

	// Clear the lattice
	void Lattice::set_all_units_empty() {
		for (auto &s: _latt) {
			s.set_b_mode_empty();
			s.set_p_mode_empty();
		};
	};

	// Binary/probabilistic
	void Lattice::convert_all_units_to_b_mode() {
		for (auto &s: _latt) {
			s.convert_p_to_b_mode();
		};
	};
	void Lattice::convert_all_units_to_p_mode() {
		for (auto &s: _latt) {
			s.convert_b_to_p_mode();
		};
	};

	/********************
	Write/read Latt to a file
	********************/

	void Lattice::write_to_file(std::string fname) 
	{
		std::ofstream f;
		f.open (fname);
		for (auto const &l: _latt) {
			if (_dim == 1) {
				f << l.x() << "\n";
			} else if (_dim == 2) {
				f << l.x() << " " << l.y() << "\n";
			} else if (_dim == 3) {
				f << l.x() << " " << l.y() << " " << l.z() << "\n";
			};
		};
		f.close();
	};
	void Lattice::read_from_file(std::string fname, bool binary)
	{
		/*
		// Clear the current Latt
		set_all_units_empty();

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
			    if (_dim == 2 || _dim == 3) {
				    iss >> y;
			    };
			    if (_dim == 3) {
				    iss >> z;
			    };
			    iss >> sp;
			    if (!binary) {
			    	// Read the prob
			    	iss >> prob;
			    };
		    	// Add to Latt
		    	if (_dim == 1) {
		    		s = _look_up(atoi(x.c_str()));
			    } else if (_dim == 2) {
		    		s = _look_up(atoi(x.c_str()),atoi(y.c_str()));
			    } else if (_dim == 3) {
			    	s = _look_up(atoi(x.c_str()),atoi(y.c_str()),atoi(z.c_str()));
			    };
			    if (binary) {
		    		s->set_prob(_sp_map[sp],1.0);
		    	} else {
		    		prob_val = atof(prob.c_str());
		    		s->set_prob(_sp_map[sp],prob_val);
		    		s->set_prob(nullptr,1.0-prob_val);
		    	};
	    		// Reset
		    	sp=""; x=""; y=""; z=""; prob="";
			};
		};
		f.close();

		// std::cout << "Read: " << fname << std::endl;
		*/
	};

	/********************
	Sample
	********************/

	void Lattice::sample_at_timepoint(int timepoint, bool binary) {

		for (auto &s: _latt) {
			s.sample_at_timepoint(timepoint, binary);
		};
	};

	/********************
	Get counts
	********************/

	// 1 particle
	double Lattice::get_count(Sptr &sp, bool binary) const {
		double count = 0.0;
		if (binary) {
			for (auto const &s: _latt) {
				if (sp == s.get_b_mode_species()) {
					count += 1.0;
				};
			};
		} else {
			for (auto const &s: _latt) {
				count += s.get_p_mode_prob(sp);
			};
		};

		return count;
	};

	// 2 particle
	void Lattice::_get_count(double &count, Sptr &sp1, Sptr &sp2, const UnitVisible &uv1, const UnitVisible *uv2, bool binary, bool reversibly) const {
		if (binary) {
			if ((sp1 == uv1.get_b_mode_species()) && (sp2 == uv2->get_b_mode_species())) {
				count += 1.0;
			};
			if (reversibly && sp1 != sp2) {
				if ((sp1 == uv2->get_b_mode_species()) && (sp2 == uv1.get_b_mode_species())) {
					count += 1.0;
				};
			};
		} else {
			count += uv1.get_p_mode_prob(sp1) * uv2->get_p_mode_prob(sp2); 
			if (reversibly && sp1 != sp2) {
				count += uv1.get_p_mode_prob(sp2) * uv2->get_p_mode_prob(sp1); 
			};
		};
	};
	double Lattice::get_count(Sptr &sp1, Sptr &sp2, bool binary, bool reversibly) const {
		exit(EXIT_FAILURE);

		const UnitVisible *nbr = nullptr;
		double count = 0.0;
		for (auto &s: _latt) {
			// Only connect to "plus one" (not minus one) => no duplicates!

			// Dim 1,2,3
			if (s.x()+1 <= _box_length) {
				if (_dim == 1) {
					nbr = &_look_up_const(s.x()+1);
				} else if (_dim == 2) {
					nbr = &_look_up_const(s.x()+1,s.y());
				} else if (_dim == 3) {
					nbr = &_look_up_const(s.x()+1,s.y(),s.z());
				};

				_get_count(count, sp1, sp2, s, nbr, binary, reversibly);
			};

			// Dim 2,3
			if (s.y()+1 <= _box_length) {
				if (_dim == 2) {
					nbr = &_look_up_const(s.x(),s.y()+1);
				} else if (_dim == 3) {
					nbr = &_look_up_const(s.x(),s.y()+1,s.z());
				};

				_get_count(count, sp1, sp2, s, nbr, binary, reversibly);
			};

			// Dim 3
			if (s.z()+1 <= _box_length) {
				if (_dim == 3) {
					nbr = &_look_up_const(s.x(),s.y(),s.z()+1);
				};

				_get_count(count, sp1, sp2, s, nbr, binary, reversibly);
			};
		};

		return count;
	};

	// 3 particle
	void Lattice::_get_count(double &count, Sptr &sp1, Sptr &sp2, Sptr &sp3, const UnitVisible &uv1, const UnitVisible *uv2, const UnitVisible *uv3, bool binary, bool reversibly) const {
		if (binary) {
			if ((sp1 == uv1.get_b_mode_species()) && (sp2 == uv2->get_b_mode_species()) && (sp3 == uv3->get_b_mode_species())) {
				count += 1.0;
			};
			if (reversibly && sp1 != sp3) {
				if ((sp3 == uv1.get_b_mode_species()) && (sp2 == uv2->get_b_mode_species()) && (sp1 == uv3->get_b_mode_species())) {
					count += 1.0;
				};
			};
		} else {
			count += uv1.get_p_mode_prob(sp1) * uv2->get_p_mode_prob(sp2) * uv3->get_p_mode_prob(sp3); 
			if (reversibly && sp1 != sp3) {
				count += uv1.get_p_mode_prob(sp3) * uv2->get_p_mode_prob(sp2) * uv3->get_p_mode_prob(sp1); 
			};
		};
	};

	double Lattice::get_count(Sptr &sp1, Sptr &sp2, Sptr &sp3, bool binary, bool reversibly) const {
		exit(EXIT_FAILURE);
		
		double count = 0.0;
		
		// ugh...

		return count;	
	};

	/****************************************
	Lattice - PRIVATE
	****************************************/

	/********************
	Lookup a site iterator from x,y,z
	********************/

	const UnitVisible& Lattice::_look_up_const(int x) const {
		// Figure out index in list
		int n = x-1;

		// Grab
		Latt::const_iterator it = _latt.begin();
		std::advance(it,n);
		return *it;
	};
	UnitVisible& Lattice::_look_up(int x) {
		// Figure out index in list
		int n = x-1;

		// Grab
		Latt_it it = _latt.begin();
		std::advance(it,n);
		return *it;
	};
	const UnitVisible& Lattice::_look_up_const(int x, int y) const {
		// Figure out index in list
		int n = (x-1)*_box_length + y-1;

		// Grab
		Latt::const_iterator it = _latt.begin();
		std::advance(it,n);
		return *it;
	};
	UnitVisible& Lattice::_look_up(int x, int y) {
		// Figure out index in list
		int n = (x-1)*_box_length + y-1;

		// Grab
		Latt_it it = _latt.begin();
		std::advance(it,n);
		return *it;
	};
	const UnitVisible& Lattice::_look_up_const(int x, int y, int z) const {
		// Figure out index in list
		int n = (x-1)*_box_length*_box_length + (y-1)*_box_length + z-1;

		// Grab
		Latt::const_iterator it = _latt.begin();
		std::advance(it,n);
		return *it;
	};
	UnitVisible& Lattice::_look_up(int x, int y, int z) {
		// Figure out index in list
		int n = (x-1)*_box_length*_box_length + (y-1)*_box_length + z-1;

		// Grab
		Latt_it it = _latt.begin();
		std::advance(it,n);
		return *it;
	};

	/********************
	Check dim
	********************/

	void Lattice::_check_dim(int dim) const {
		if (_dim != dim) {
			std::cerr << ">>> Error: Lattice::_check_dim <<< dim is: " << _dim << " but requested: " << dim << std::endl;
			exit(EXIT_FAILURE);
		};	
	};
};