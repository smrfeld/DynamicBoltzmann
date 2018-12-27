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

	double ConnVV::get_moment(std::string ixn_param_name, bool binary) const {
		if (!_ixn_dict) {
			std::cerr << ">>> Error: ConnVV::get_moment <<< no ixn dict exists on this conn" << std::endl;
			exit(EXIT_FAILURE);
		};

		double count = 0.0;

		std::vector<Sptr2> species = _ixn_dict->get_species_from_ixn(ixn_param_name);
		for (auto &sp: species) {

			if (binary) {

				// Binary
				if (_uv1->get_b_mode_species() == sp.s1 && _uv2->get_b_mode_species() == sp.s2) {
					count += 1.0;
				};

			} else {

				// Prob
				count += _uv1->get_p_mode_prob(sp.s1) * _uv2->get_p_mode_prob(sp.s2);
			};
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
		if (v_other->check_is_b_mode()) {
			//////////
			// Binary
			/////////
			Sptr sp_other = v_other->get_b_mode_species();
			// Get ixn
			if (idx==0) {
				ixn = _ixn_dict->get_ixn(sp_to_place,sp_other);
			} else if (idx==1) {
				ixn = _ixn_dict->get_ixn(sp_other,sp_to_place);
			};
			// Add
			act += ixn; // * 1.0 * 1.0
		} else {
			////////////////
			// Probabilistic
			////////////////
			for (auto const &pr: v_other->get_p_mode_probs()) {
				if (pr.second == 0.0) {
					continue;
				};
				// Get ixn
				if (idx==0) {
					ixn = _ixn_dict->get_ixn(sp_to_place,pr.first);
				} else if (idx==1) {
					ixn = _ixn_dict->get_ixn(pr.first,sp_to_place);
				};
				// Add
				act += ixn * pr.second; // * 1.0
			};
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
		// Nothing....
	};
	void ConnVH::_move(ConnVH& other) {
		_uv = std::move(other._uv);
		_uh = std::move(other._uh);
		_ixn_dict = std::move(other._ixn_dict);
	};
	void ConnVH::_copy(const ConnVH& other) {
		_uv = other._uv;
		_uh = other._uh;
		_ixn_dict = other._ixn_dict;
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

	double ConnVH::get_moment(std::string ixn_param_name, bool binary_visible, bool binary_hidden) const {
		if (!_ixn_dict) {
			std::cerr << ">>> Error: ConnVH::get_moment <<< no ixn dict exists on this conn" << std::endl;
			exit(EXIT_FAILURE);
		};

		if (binary_visible && binary_hidden) {
			// Binary - binary
			int count = 0;
			for (auto const &sp: _ixn_dict->get_species_from_ixn(ixn_param_name)) {
				if (_uv->check_is_b_mode_species(sp.s1) && _uh->check_is_b_mode_species(sp.s2)) {
					count++;
				};
			};
			return count;
		} else if (binary_visible && !binary_hidden) {
			// Binary - prob
			double count = 0.0;
			for (auto const &sp: _ixn_dict->get_species_from_ixn(ixn_param_name)) {
				if (_uv->check_is_b_mode_species(sp.s1)) {
					count += 1.0 * _uh->get_p_mode_prob(sp.s2);
				};
			};
			return count;
		} else if (!binary_visible && binary_hidden) {
			// Prob - binary
			double count = 0.0;
			for (auto const &sp: _ixn_dict->get_species_from_ixn(ixn_param_name)) {
				if (_uh->check_is_b_mode_species(sp.s2)) {
					count += 1.0 * _uv->get_p_mode_prob(sp.s1);
				};
			};
			return count;
		} else {
			// Prob - prob
			double count = 0.0;
			for (auto const &sp: _ixn_dict->get_species_from_ixn(ixn_param_name)) {
				count += _uv->get_p_mode_prob(sp.s1) * _uh->get_p_mode_prob(sp.s2);
			};
			return count;
		};
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
		if (_uh->check_is_b_mode()) {
			//////////
			// Binary
			/////////
			Sptr sp_other = _uh->get_b_mode_species();
			// Get ixn
			ixn = _ixn_dict->get_ixn(sp_to_place,sp_other);
			// Add
			act += ixn; // * 1.0 * 1.0
		} else {
			////////////////
			// Probabilistic
			////////////////
			for (auto const &pr: _uh->get_p_mode_probs()) {
				if (pr.second == 0.0) {
					continue;
				};
				// Get ixn
				ixn = _ixn_dict->get_ixn(sp_to_place,pr.first);
				// Add
				act += ixn * pr.second; // * 1.0
			};
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
		if (_uv->check_is_b_mode()) {
			//////////
			// Binary
			/////////
			Sptr sp_other = _uv->get_b_mode_species();
			// Get ixn
			ixn = _ixn_dict->get_ixn(sp_other,sp_to_place);
			// Add
			act += ixn; // * 1.0 * 1.0
		} else {
			////////////////
			// Probabilistic
			////////////////
			for (auto const &pr: _uv->get_p_mode_probs()) {
				if (pr.second == 0.0) {
					continue;
				};
				// Get ixn
				ixn = _ixn_dict->get_ixn(pr.first,sp_to_place);
				// Add
				act += ixn * pr.second; // * 1.0
			};
		};

		return act;
	};





























	/****************************************
	Class to hold a connection
	****************************************/

	// Constructor
	ConnHH::ConnHH(UnitHidden *uh1, UnitHidden *uh2) {
		_uh1 = uh1;
		_uh2 = uh2;
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
		_uh1 = std::move(other._uh1);
		_uh2 = std::move(other._uh2);
		_ixn_dict = std::move(other._ixn_dict);
	};
	void ConnHH::_copy(const ConnHH& other) {
		_uh1 = other._uh1;
		_uh2 = other._uh2;
		_ixn_dict = other._ixn_dict;
	};

	/********************
	Getters
	********************/

	UnitHidden* ConnHH::get_unit_h(int idx) const {
		if (idx == 0) {
			return _uh1;
		} else {
			return _uh2;
		};
	};

	/********************
	Check units involved
	********************/

	bool ConnHH::check_connects_unit(UnitHidden *uh) const {
		if (_uh1 == uh || _uh2 == uh) {
			return true;
		};
		return false;
	};
	bool ConnHH::check_connects_units(UnitHidden *uh1, UnitHidden *uh2, bool this_order) const {
		if (this_order) {
			if (_uh1 == uh1 && _uh2 == uh2) {
				return true;
			};
		} else {
			if ((_uh1 == uh1 && _uh2 == uh2) || (_uh1 == uh2 && _uh2 == uh1)) {
				return true;
			};
		};
		return false;
	};

	/********************
	Get count
	********************/

	double ConnHH::get_moment(std::string ixn_param_name, bool binary) const {
		if (!_ixn_dict) {
			std::cerr << ">>> Error: ConnHH::get_moment <<< no ixn dict exists on this conn" << std::endl;
			exit(EXIT_FAILURE);
		};

		if (binary) {
			// Binary
			int count = 0;
			for (auto const &sp: _ixn_dict->get_species_from_ixn(ixn_param_name)) {
				if (_uh1->check_is_b_mode_species(sp.s1) && _uh2->check_is_b_mode_species(sp.s2)) {
					count++;
				};
			};
			return count;
		} else {
			// Prob
			double count = 0.0;
			for (auto const &sp: _ixn_dict->get_species_from_ixn(ixn_param_name)) {
				count += _uh1->get_p_mode_prob(sp.s1) * _uh2->get_p_mode_prob(sp.s2);
			};
			return count;
		};
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
			h_other = _uh2;
		} else if (idx==1) {
			h_other = _uh1;
		};

		double ixn;
		// Go through species of other site
		if (h_other->check_is_b_mode()) {
			//////////
			// Binary
			/////////
			Sptr sp_other = h_other->get_b_mode_species();
			// Get ixn
			if (idx==0) {
				ixn = _ixn_dict->get_ixn(sp_to_place,sp_other);
			} else if (idx==1) {
				ixn = _ixn_dict->get_ixn(sp_other,sp_to_place);
			};
			// Add
			act += ixn; // * 1.0 * 1.0
		} else {
			////////////////
			// Probabilistic
			////////////////
			for (auto const &pr: h_other->get_p_mode_probs()) {
				if (pr.second == 0.0) {
					continue;
				};
				// Get ixn
				if (idx==0) {
					ixn = _ixn_dict->get_ixn(sp_to_place,pr.first);
				} else if (idx==1) {
					ixn = _ixn_dict->get_ixn(pr.first,sp_to_place);
				};
				// Add
				act += ixn * pr.second; // * 1.0
			};
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

	double ConnVVV::get_moment(std::string ixn_param_name, bool binary) const {
		if (!_ixn_dict) {
			std::cerr << ">>> Error: ConnVVV::get_moment <<< no ixn dict exist on this conn" << std::endl;
			exit(EXIT_FAILURE);
		};

		double count = 0.0;

		std::vector<Sptr3> species = _ixn_dict->get_species_from_ixn(ixn_param_name);
		for (auto &sp: species) {

			if (binary) {

				// Binary
				if (_uv1->get_b_mode_species() == sp.s1 && _uv2->get_b_mode_species() == sp.s2 && _uv3->get_b_mode_species() == sp.s3) {
					count += 1.0;
				};

			} else {

				// Prob
				count += _uv1->get_p_mode_prob(sp.s1) * _uv2->get_p_mode_prob(sp.s2) * _uv3->get_p_mode_prob(sp.s3);
			};
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
		if (v_other_1->check_is_b_mode()) {

			//////////
			// Binary
			/////////

			Sptr sp_other_1 = v_other_1->get_b_mode_species();

			// Go through species of other site #2
			if (v_other_2->check_is_b_mode()) {

				//////////
				// Binary - Binary
				/////////

				Sptr sp_other_2 = v_other_2->get_b_mode_species();

				// Get ixn
				if (idx==0) {
					ixn = _ixn_dict->get_ixn(sp_to_place,sp_other_1,sp_other_2);
				} else if (idx==1) {
					ixn = _ixn_dict->get_ixn(sp_other_1,sp_to_place,sp_other_2);
				} else if (idx==2) {
					ixn = _ixn_dict->get_ixn(sp_other_1,sp_other_2,sp_to_place);
				};

				// Add
				act += ixn;

			} else {

				//////////
				// Binary - Probabilistic
				/////////

				for (auto const &pr: v_other_2->get_p_mode_probs()) {
					if (pr.second == 0.0) {
						continue;
					};
					// Get ixn
					if (idx==0) {
						ixn = _ixn_dict->get_ixn(sp_to_place,sp_other_1,pr.first);
					} else if (idx==1) {
						ixn = _ixn_dict->get_ixn(sp_other_1,sp_to_place,pr.first);
					} else if (idx==2) {
						ixn = _ixn_dict->get_ixn(sp_other_1,pr.first,sp_to_place);
					};

					// Add
					act += ixn * pr.second;
				};

			};

		} else {

			////////////////
			// Probabilistic
			////////////////

			for (auto const &pr: v_other_1->get_p_mode_probs()) {
				if (pr.second == 0.0) {
					continue;
				};

				// Go through species of other site #2
				if (v_other_2->check_is_b_mode()) {

					//////////
					// Probabilistic - Binary
					/////////

					Sptr sp_other_2 = v_other_2->get_b_mode_species();

					// Get ixn
					if (idx==0) {
						ixn = _ixn_dict->get_ixn(sp_to_place,pr.first,sp_other_2);
					} else if (idx==1) {
						ixn = _ixn_dict->get_ixn(pr.first,sp_to_place,sp_other_2);
					} else if (idx==2) {
						ixn = _ixn_dict->get_ixn(pr.first,sp_other_2,sp_to_place);
					};
				
					// Add
					act += ixn * pr.second;

				} else {

					//////////
					// Probabilistic - Probabilistic
					/////////

					for (auto const &pr2: v_other_2->get_p_mode_probs()) {
						if (pr2.second == 0.0) {
							continue;
						};
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
			};
		};

		return act;
	};

};