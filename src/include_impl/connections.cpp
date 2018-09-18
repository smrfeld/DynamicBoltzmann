#include "../../include/dynamicboltz_bits/connections.hpp"

// Other headers
#include "../../include/dynamicboltz_bits/general.hpp"
#include "../../include/dynamicboltz_bits/unit_visible.hpp"
#include "../../include/dynamicboltz_bits/ixn_dicts.hpp"

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

	double ConnVV::get_count(Sptr &sp1, Sptr &sp2, bool binary, bool this_order) const {
		double count = 0.0;

		if (binary) {

			// Binary

			if (_uv1->get_b_mode_species() == sp1 && _uv2->get_b_mode_species() == sp2) {
				count += 1.0;
			};
			if (!this_order && (sp1 != sp2)) {
				if (_uv1->get_b_mode_species() == sp2 && _uv2->get_b_mode_species() == sp1) {
					count += 1.0;
				};
			};

		} else {

			// Prob

			count += _uv1->get_p_mode_prob(sp1) * _uv2->get_p_mode_prob(sp2);
			if (!this_order && (sp1 != sp2)) {
				count += _uv1->get_p_mode_prob(sp2) * _uv2->get_p_mode_prob(sp1);
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

	double ConnVV::get_act_for_species_at_unit_at_timepoint(Sptr &sp_to_place, int idx, int timepoint) {
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
				ixn = _ixn_dict->get_ixn_at_timepoint(sp_to_place,sp_other,timepoint);
			} else if (idx==1) {
				ixn = _ixn_dict->get_ixn_at_timepoint(sp_other,sp_to_place,timepoint);
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
					ixn = _ixn_dict->get_ixn_at_timepoint(sp_to_place,pr.first,timepoint);
				} else if (idx==1) {
					ixn = _ixn_dict->get_ixn_at_timepoint(pr.first,sp_to_place,timepoint);
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

	double ConnVVV::get_count(Sptr &sp1, Sptr &sp2, Sptr &sp3, bool binary, bool this_order) const {
		double count = 0.0;

		if (binary) {

			// Binary

			if (_uv1->get_b_mode_species() == sp1 && _uv2->get_b_mode_species() == sp2 && _uv3->get_b_mode_species() == sp3) {
				count += 1.0;
			};
			if (!this_order && sp1 != sp3) {
				if (_uv1->get_b_mode_species() == sp3 && _uv2->get_b_mode_species() == sp2 && _uv3->get_b_mode_species() == sp1) {
					count += 1.0;
				};
			};

		} else {

			// Prob

			count += _uv1->get_p_mode_prob(sp1) * _uv2->get_p_mode_prob(sp2) * _uv3->get_p_mode_prob(sp3);
			if (!this_order && (sp1 != sp3)) {
				count += _uv1->get_p_mode_prob(sp3) * _uv2->get_p_mode_prob(sp2) * _uv3->get_p_mode_prob(sp1);
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

	double ConnVVV::get_act_for_species_at_unit_at_timepoint(Sptr &sp_to_place, int idx, int timepoint) {
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
					ixn = _ixn_dict->get_ixn_at_timepoint(sp_to_place,sp_other_1,sp_other_2,timepoint);
				} else if (idx==1) {
					ixn = _ixn_dict->get_ixn_at_timepoint(sp_other_1,sp_to_place,sp_other_2,timepoint);
				} else if (idx==2) {
					ixn = _ixn_dict->get_ixn_at_timepoint(sp_other_1,sp_other_2,sp_to_place,timepoint);
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
						ixn = _ixn_dict->get_ixn_at_timepoint(sp_to_place,sp_other_1,pr.first,timepoint);
					} else if (idx==1) {
						ixn = _ixn_dict->get_ixn_at_timepoint(sp_other_1,sp_to_place,pr.first,timepoint);
					} else if (idx==2) {
						ixn = _ixn_dict->get_ixn_at_timepoint(sp_other_1,pr.first,sp_to_place,timepoint);
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
						ixn = _ixn_dict->get_ixn_at_timepoint(sp_to_place,pr.first,sp_other_2,timepoint);
					} else if (idx==1) {
						ixn = _ixn_dict->get_ixn_at_timepoint(pr.first,sp_to_place,sp_other_2,timepoint);
					} else if (idx==2) {
						ixn = _ixn_dict->get_ixn_at_timepoint(pr.first,sp_other_2,sp_to_place,timepoint);
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
							ixn = _ixn_dict->get_ixn_at_timepoint(sp_to_place,pr.first,pr2.first,timepoint);
						} else if (idx==1) {
							ixn = _ixn_dict->get_ixn_at_timepoint(pr.first,sp_to_place,pr2.first,timepoint);
						} else if (idx==2) {
							ixn = _ixn_dict->get_ixn_at_timepoint(pr.first,pr2.first,sp_to_place,timepoint);
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