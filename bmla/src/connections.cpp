#include "../include/bmla_bits/connections.hpp"

// Other headers
#include "../include/bmla_bits/general.hpp"
#include "../include/bmla_bits/unit.hpp"
#include "../include/bmla_bits/ixn_dicts.hpp"
#include "../include/bmla_bits/species.hpp"

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
	Class to hold a connection
	****************************************/

	// Constructor
	ConnVV::ConnVV(UnitVisible *uv1, UnitVisible *uv2) {
		_uv1 = uv1;
		_uv2 = uv2;
	};
	ConnVV::ConnVV(const ConnVV& other) {
		_copy(other);
	};
	ConnVV::ConnVV(ConnVV&& other) {
		_move(other);
	};
	ConnVV& ConnVV::operator=(const ConnVV& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	ConnVV& ConnVV::operator=(ConnVV&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
	};
	ConnVV::~ConnVV() {
		_clean_up();
	};
	void ConnVV::_clean_up() {
		// Nothing....
	};
	void ConnVV::_move(ConnVV& other) {
		_uv1 = std::move(other._uv1);
		_uv2 = std::move(other._uv2);
		_ixn_dict = std::move(other._ixn_dict);
	};
	void ConnVV::_copy(const ConnVV& other) {
		_uv1 = other._uv1;
		_uv2 = other._uv2;
		_ixn_dict = other._ixn_dict;
	};

	/********************
	Check units involved
	********************/

	bool ConnVV::check_connects_unit(UnitVisible *uv) const {
		if (_uv1 == uv || _uv2 == uv) {
			return true;
		};
		return false;
	};
	bool ConnVV::check_connects_units(UnitVisible *uv1, UnitVisible *uv2, bool this_order) const {
		if (this_order) {
			if (_uv1 == uv1 && _uv2 == uv2) {
				return true;
			};
		} else {
			if ((_uv1 == uv1 && _uv2 == uv2) || (_uv1 == uv2 && _uv2 == uv1)) {
				return true;
			};
		};
		return false;
	};

	/********************
	Get count
	********************/

	double ConnVV::get_moment(std::string ixn_param_name) const {
		if (!_ixn_dict) {
			std::cerr << ">>> Error: ConnVV::get_moment <<< no ixn dict exists on this conn" << std::endl;
			exit(EXIT_FAILURE);
		};

		double count = 0.0;

		std::vector<Sptr2> species = _ixn_dict->get_species_from_ixn(ixn_param_name);
		for (auto &sp: species) {
			count += _uv1->get_occ(sp.s1) * _uv2->get_occ(sp.s2);
		};

		return count;
	};

	/********************
	Set ixns
	********************/

	void ConnVV::set_ixn_dict(std::shared_ptr<O2IxnDict> &ixn_dict) {
		_ixn_dict = ixn_dict;
	};

	/********************
	Get activation on site
	********************/

	double ConnVV::get_act_for_species_at_unit(const Sptr &sp_to_place, int idx) {
		double act=0.0;

		// 0 if no ixn dict
		if (!_ixn_dict) {
			return act;
		};

		UnitVisible *v_other=nullptr;
		if (idx==0) {
			v_other = _uv2;
		} else if (idx==1) {
			v_other = _uv1;
		};

		double ixn;
		// Go through species of other site
		for (auto const &pr: v_other->get_nonzero_occs()) {
			// Get ixn
			if (idx==0) {
				ixn = _ixn_dict->get_ixn(sp_to_place,pr.first);
			} else if (idx==1) {
				ixn = _ixn_dict->get_ixn(pr.first,sp_to_place);
			};
			// Add
			act += ixn * pr.second; // * 1.0
		};

		return act;
	};








































	/****************************************
	Class to hold a connection
	****************************************/

	// Constructor
	ConnVH::ConnVH(UnitVisible *uv, UnitHidden* uh) {
		_uv = uv;
		_uh = uh;
		_multiplier_h_to_v = nullptr;
        _multiplier_v_to_h = nullptr;
	};
	ConnVH::ConnVH(const ConnVH& other) {
		_copy(other);
	};
	ConnVH::ConnVH(ConnVH&& other) {
		_move(other);
	};
	ConnVH& ConnVH::operator=(const ConnVH& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	ConnVH& ConnVH::operator=(ConnVH&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
	};
	ConnVH::~ConnVH() {
		_clean_up();
	};
	void ConnVH::_clean_up() {
        if (_multiplier_h_to_v) {
            delete _multiplier_h_to_v;
            _multiplier_h_to_v = nullptr;
        };
		if (_multiplier_v_to_h) {
			delete _multiplier_v_to_h;
            _multiplier_v_to_h = nullptr;
		};
	};
	void ConnVH::_move(ConnVH& other) {
		_uv = std::move(other._uv);
		_uh = std::move(other._uh);
		_ixn_dict = std::move(other._ixn_dict);
		_multiplier_h_to_v = other._multiplier_h_to_v;
        _multiplier_v_to_h = other._multiplier_v_to_h;
        
		other._multiplier_h_to_v = nullptr;
        other._multiplier_v_to_h = nullptr;
	};
	void ConnVH::_copy(const ConnVH& other) {
		_uv = other._uv;
		_uh = other._uh;
		_ixn_dict = other._ixn_dict;
		if (other._multiplier_h_to_v) {
			_multiplier_h_to_v = new double(*other._multiplier_h_to_v);
        };
        if (other._multiplier_v_to_h) {
            _multiplier_v_to_h = new double(*other._multiplier_v_to_h);
		};	
	};

	/********************
	Multiplier
	********************/

	void ConnVH::set_multiplier_h_to_v(double multiplier_h_to_v) {
		if (!_multiplier_h_to_v) {
			_multiplier_h_to_v = new double(multiplier_h_to_v);
		} else {
			*_multiplier_h_to_v = multiplier_h_to_v;
		};
	};
    void ConnVH::set_multiplier_v_to_h(double multiplier_v_to_h) {
        if (!_multiplier_v_to_h) {
            _multiplier_v_to_h = new double(multiplier_v_to_h);
        } else {
            *_multiplier_v_to_h = multiplier_v_to_h;
        };
    };
    
	/********************
	Check units involved
	********************/

	bool ConnVH::check_connects_unit(UnitVisible *uv) const {
		if (_uv == uv) {
			return true;
		};
		return false;
	};
	bool ConnVH::check_connects_unit(UnitHidden *uh) const {
		if (_uh == uh) {
			return true;
		};
		return false;
	};
	bool ConnVH::check_connects_units(UnitVisible *uv, UnitHidden *uh) const {
		if (_uv == uv && _uh == uh) {
			return true;
		};
		return false;
	};

	/********************
	Get units
	********************/

	UnitVisible* ConnVH::get_unit_v() const {
		return _uv;
	};
	UnitHidden* ConnVH::get_unit_h() const {
		return _uh;
	};

	/********************
	Get count
	********************/

	double ConnVH::get_moment(std::string ixn_param_name) const {
		if (!_ixn_dict) {
			std::cerr << ">>> Error: ConnVH::get_moment <<< no ixn dict exists on this conn" << std::endl;
			exit(EXIT_FAILURE);
		};

		// Prob - prob
		double count = 0.0;
		for (auto const &sp: _ixn_dict->get_species_from_ixn(ixn_param_name)) {
			count += _uv->get_occ(sp.s1) * _uh->get_occ(sp.s2);
		};
		return count;
	};

	/********************
	Set ixns
	********************/

	void ConnVH::set_ixn_dict(std::shared_ptr<O2IxnDict> &ixn_dict) {
		_ixn_dict = ixn_dict;
	};

	/********************
	Get activation on site
	********************/

	double ConnVH::get_act_for_species_at_unit_v(const Sptr &sp_to_place) {
		double act=0.0;

		// 0 if no ixn dict
		if (!_ixn_dict) {
			return act;
		};

		double ixn;
		// Go through species of other site
		for (auto const &pr: _uh->get_nonzero_occs()) {
			// Get ixn
			ixn = _ixn_dict->get_ixn(sp_to_place,pr.first);
			// Add
			act += ixn * pr.second; // * 1.0
		};

		if (_multiplier_h_to_v) {
			act *= *_multiplier_h_to_v;
		};

		return act;
	};

	double ConnVH::get_act_for_species_at_unit_h(const Sptr &sp_to_place) {
		double act=0.0;

		// 0 if no ixn dict
		if (!_ixn_dict) {
			return act;
		};

		double ixn;
		// Go through species of other site
		for (auto const &pr: _uv->get_nonzero_occs()) {
			// Get ixn
			ixn = _ixn_dict->get_ixn(pr.first,sp_to_place);
			// Add
			act += ixn * pr.second; // * 1.0
		};

		if (_multiplier_v_to_h) {
			act *= *_multiplier_v_to_h;
		};

		return act;
	};





























	/****************************************
	Class to hold a connection
	****************************************/

	// Constructor
	ConnHH::ConnHH(UnitHidden *uh1, UnitHidden *uh2) {
		_uh0 = uh1;
		_uh1 = uh2;
        _multiplier_idx_0 = nullptr;
        _multiplier_idx_1 = nullptr;
    };
	ConnHH::ConnHH(const ConnHH& other) {
		_copy(other);
	};
	ConnHH::ConnHH(ConnHH&& other) {
		_move(other);
	};
	ConnHH& ConnHH::operator=(const ConnHH& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	ConnHH& ConnHH::operator=(ConnHH&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
	};
	ConnHH::~ConnHH() {
		_clean_up();
	};
	void ConnHH::_clean_up() {
		// Nothing....
	};
	void ConnHH::_move(ConnHH& other) {
		_uh0 = std::move(other._uh0);
		_uh1 = std::move(other._uh1);
		_ixn_dict = std::move(other._ixn_dict);
        _multiplier_idx_0 = other._multiplier_idx_0;
        _multiplier_idx_1 = other._multiplier_idx_1;

        other._multiplier_idx_0 = nullptr;
        other._multiplier_idx_1 = nullptr;
	};
	void ConnHH::_copy(const ConnHH& other) {
		_uh0 = other._uh0;
		_uh1 = other._uh1;
		_ixn_dict = other._ixn_dict;
        if (other._multiplier_idx_0) {
            _multiplier_idx_0 = new double(*other._multiplier_idx_0);
        };
        if (other._multiplier_idx_1) {
            _multiplier_idx_1 = new double(*other._multiplier_idx_1);
        };
	};

    /********************
     Multiplier
     ********************/
    
    void ConnHH::set_multiplier_for_unit(int idx, double multiplier) {
        if (idx == 0) {
            if (!_multiplier_idx_0) {
                _multiplier_idx_0 = new double(multiplier);
            } else {
                *_multiplier_idx_0 = multiplier;
            };
        } else {
            if (!_multiplier_idx_1) {
                _multiplier_idx_1 = new double(multiplier);
            } else {
                *_multiplier_idx_1 = multiplier;
            };
        };
    };
    
	/********************
	Getters
	********************/

	UnitHidden* ConnHH::get_unit_h(int idx) const {
		if (idx == 0) {
			return _uh0;
		} else {
			return _uh1;
		};
	};

	/********************
	Check units involved
	********************/

	bool ConnHH::check_connects_unit(UnitHidden *uh) const {
		if (_uh0 == uh || _uh1 == uh) {
			return true;
		};
		return false;
	};
	bool ConnHH::check_connects_units(UnitHidden *uh1, UnitHidden *uh2, bool this_order) const {
		if (this_order) {
			if (_uh0 == uh1 && _uh1 == uh2) {
				return true;
			};
		} else {
			if ((_uh0 == uh1 && _uh1 == uh2) || (_uh0 == uh2 && _uh1 == uh1)) {
				return true;
			};
		};
		return false;
	};

	/********************
	Get count
	********************/

	double ConnHH::get_moment(std::string ixn_param_name) const {
		if (!_ixn_dict) {
			std::cerr << ">>> Error: ConnHH::get_moment <<< no ixn dict exists on this conn" << std::endl;
			exit(EXIT_FAILURE);
		};

		// Prob
		double count = 0.0;
		for (auto const &sp: _ixn_dict->get_species_from_ixn(ixn_param_name)) {
			count += _uh0->get_occ(sp.s1) * _uh1->get_occ(sp.s2);
		};
		return count;
	};

	/********************
	Set ixns
	********************/

	void ConnHH::set_ixn_dict(std::shared_ptr<O2IxnDict> &ixn_dict) {
		_ixn_dict = ixn_dict;
	};

	/********************
	Get activation on site
	********************/

	double ConnHH::get_act_for_species_at_unit(const Sptr &sp_to_place, int idx) {
		double act=0.0;

		// 0 if no ixn dict
		if (!_ixn_dict) {
			return act;
		};

		UnitHidden *h_other=nullptr;
		if (idx==0) {
			h_other = _uh1;
		} else if (idx==1) {
			h_other = _uh0;
		};

		double ixn;
		// Go through species of other site
		for (auto const &pr: h_other->get_nonzero_occs()) {
			// Get ixn
			if (idx==0) {
				ixn = _ixn_dict->get_ixn(sp_to_place,pr.first);
			} else if (idx==1) {
				ixn = _ixn_dict->get_ixn(pr.first,sp_to_place);
			};
			// Add
			act += ixn * pr.second; // * 1.0
		};
        
        if (idx == 0 && _multiplier_idx_0) {
            act *= *_multiplier_idx_0;
        } else if (idx == 1 && _multiplier_idx_1) {
            act *= *_multiplier_idx_1;
        };

		return act;
	};











































	/****************************************
	Class to hold a connection
	****************************************/

	// Constructor
	ConnVVV::ConnVVV(UnitVisible *uv1, UnitVisible *uv2, UnitVisible *uv3) {
		_uv1 = uv1;
		_uv2 = uv2;
		_uv3 = uv3;
	};
	ConnVVV::ConnVVV(const ConnVVV& other) {
		_copy(other);
	};
	ConnVVV::ConnVVV(ConnVVV&& other) {
		_move(other);
	};
	ConnVVV& ConnVVV::operator=(const ConnVVV& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	ConnVVV& ConnVVV::operator=(ConnVVV&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
	};
	ConnVVV::~ConnVVV() {
		_clean_up();
	};
	void ConnVVV::_clean_up() {
		// Nothing....
	};
	void ConnVVV::_move(ConnVVV& other) {
		_uv1 = std::move(other._uv1);
		_uv2 = std::move(other._uv2);
		_uv3 = std::move(other._uv3);
		_ixn_dict = std::move(other._ixn_dict);
	};
	void ConnVVV::_copy(const ConnVVV& other) {
		_uv1 = other._uv1;
		_uv2 = other._uv2;
		_uv3 = other._uv3;
		_ixn_dict = other._ixn_dict;
	};

	/********************
	Check units involved
	********************/

	bool ConnVVV::check_connects_unit(UnitVisible *uv) {
		if (_uv1 == uv || _uv2 == uv || _uv3 == uv) {
			return true;
		};
		return false;
	};
	bool ConnVVV::check_connects_units(UnitVisible *uv1, UnitVisible *uv2, bool this_order) {
		if (this_order) {
			if ((_uv1 == uv1 && _uv2 == uv2) || (_uv2 == uv1 && _uv3 == uv2)) {
				return true;
			};
		} else {
			if ((_uv1 == uv1 && _uv2 == uv2) || (_uv1 == uv2 && _uv2 == uv1) || (_uv2 == uv1 && _uv3 == uv2) || (_uv2 == uv2 && _uv3 == uv1) || (_uv1 == uv1 && _uv3 == uv2) || (_uv1 == uv2 && _uv3 == uv1)) {
				return true;
			};
		};
		return false;
	};
	bool ConnVVV::check_connects_units(UnitVisible *uv1, UnitVisible *uv2, UnitVisible *uv3, bool this_order) {
		if (this_order) {
			if (_uv1 == uv1 && _uv2 == uv2 && _uv3 == uv3) {
				return true;
			};
		} else {
			if ((_uv1 == uv1 && _uv2 == uv2 && _uv3 == uv3) || (_uv1 == uv1 && _uv2 == uv3 && _uv3 == uv2) || (_uv1 == uv2 && _uv2 == uv1 && _uv3 == uv3) || (_uv1 == uv2 && _uv2 == uv3 && _uv3 == uv1) || (_uv1 == uv3 && _uv2 == uv1 && _uv3 == uv2) || (_uv1 == uv3 && _uv2 == uv2 && _uv3 == uv1)) {
				return true;
			};
		};
		return false;
	};

	/********************
	Get count
	********************/

	double ConnVVV::get_moment(std::string ixn_param_name) const {
		if (!_ixn_dict) {
			std::cerr << ">>> Error: ConnVVV::get_moment <<< no ixn dict exist on this conn" << std::endl;
			exit(EXIT_FAILURE);
		};

		double count = 0.0;

		std::vector<Sptr3> species = _ixn_dict->get_species_from_ixn(ixn_param_name);

		for (auto &sp: species) {
			count += _uv1->get_occ(sp.s1) * _uv2->get_occ(sp.s2) * _uv3->get_occ(sp.s3);
		};

		return count;
	};

	/********************
	Set ixns
	********************/

	void ConnVVV::set_ixn_dict(std::shared_ptr<O3IxnDict> &ixn_dict) {
		_ixn_dict = ixn_dict;
	};

	/********************
	Get activation on site
	********************/

	double ConnVVV::get_act_for_species_at_unit(const Sptr &sp_to_place, int idx) {
		double act=0.0;

		// 0 if no ixn dict
		if (!_ixn_dict) {
			return act;
		};

		UnitVisible *v_other_1=nullptr, *v_other_2 = nullptr;
		if (idx==0) {
			v_other_1 = _uv2;
			v_other_2 = _uv3;
		} else if (idx==1) {
			v_other_1 = _uv1;
			v_other_2 = _uv3;
		} else if (idx==2) {
			v_other_1 = _uv1;
			v_other_2 = _uv2;
		};

		double ixn;
		// Go through species of other site #1
		for (auto const &pr: v_other_1->get_nonzero_occs()) {
			for (auto const &pr2: v_other_2->get_nonzero_occs()) {
				// Get ixn
				if (idx==0) {
					ixn = _ixn_dict->get_ixn(sp_to_place,pr.first,pr2.first);
				} else if (idx==1) {
					ixn = _ixn_dict->get_ixn(pr.first,sp_to_place,pr2.first);
				} else if (idx==2) {
					ixn = _ixn_dict->get_ixn(pr.first,pr2.first,sp_to_place);
				};

				// Add
				act += ixn * pr.second * pr2.second;
			};
		};

		return act;
	};

};
