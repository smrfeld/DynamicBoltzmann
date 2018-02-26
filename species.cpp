#include "species.hpp"
#include <iostream>

/************************************
* Namespace for Gillespie3D
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Species
	****************************************/
	
	Species::Species(std::string nameIn) {
		name = nameIn;
		count = 0;

		// Ptrs
		_soln_traj_ptr = nullptr;
		_t_opt_ptr = nullptr;
		_h_index = -1;
	};

	/********************
	Accessor h,j
	********************/

	double Species::h() {
		if (_h_index != -1) {
			return (*_soln_traj_ptr)[*_t_opt_ptr][_h_index];
		} else {
			return 0.0;
		};
	};
	double Species::j(Species *other) {
		if (_j_index[other] != -1) {
			return (*_soln_traj_ptr)[*_t_opt_ptr][_j_index[other]];
		} else {
			return 0.0;
		};
	};

	/********************
	Comparator
	********************/

	bool operator <(const Species& a, const Species& b) {
		return a.name < b.name;
	};

};