#include "../include/dynamicboltz_bits/diff_eq_rhs.hpp"

// Other headers
#include "../include/dynamicboltz_bits/general.hpp"
#include "../include/dynamicboltz_bits/ixn_param.hpp"
#include "../include/dynamicboltz_bits/adjoint.hpp"

#include <iostream>
#include "math.h"
#include <fstream>
#include <iomanip>

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	Domain1D
	****************************************/

	Domain1D::Domain1D(Iptr ixn_param, double min, double max, int no_pts) : Dimension1D(min,max,no_pts) {
		_ixn_param = ixn_param;
	};
	Domain1D::Domain1D(const Domain1D& other) : Dimension1D(other) {
		_copy(other);
	};
	Domain1D::Domain1D(Domain1D&& other) : Dimension1D(std::move(other)) {
		_copy(other);
		other._reset();
	};
	Domain1D& Domain1D::operator=(const Domain1D& other) {
		if (this != &other)
		{
	        Dimension1D::operator=(other);
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	Domain1D& Domain1D::operator=(Domain1D&& other) {
		if (this != &other)
		{
	        Dimension1D::operator=(std::move(other));
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
	};
	Domain1D::~Domain1D() {
		_clean_up();
	};
	void Domain1D::_copy(const Domain1D& other)
	{
		_ixn_param = other._ixn_param;
	};
	void Domain1D::_reset()
	{
		_ixn_param = nullptr;
	};
	void Domain1D::_clean_up() {
	};

	/********************
	Getters
	********************/

	std::string Domain1D::get_name() const {
		return _ixn_param->get_name();
	};
	Iptr Domain1D::get_ixn_param() const {
		return _ixn_param;
	};












































	/****************************************
	DiffEqRHS
	****************************************/

	/********************
	Constructor
	********************/

	DiffEqRHS::DiffEqRHS(std::string name, Iptr parent_ixn_param, std::vector<Domain1D*> domain) : q3c1::Grid(std::vector<q3c1::Dimension1D*>(domain.begin(),domain.end())) {
		_name = name;
		_domain = domain;
		_no_dims = _domain.size();
		_parent_ixn_param = parent_ixn_param;

		if (_no_dims == 0 || _no_dims > 3) {
			std::cerr << ">>> Error: DiffEqRHS::DiffEqRHS <<< Only dims 1,2,3 are supported" << std::endl;
			exit(EXIT_FAILURE);
		};
	};
	DiffEqRHS::DiffEqRHS(const DiffEqRHS& other) : Grid(other) {
		_copy(other);
	};
	DiffEqRHS::DiffEqRHS(DiffEqRHS&& other) : Grid(std::move(other)) {
		_move(other);
	};

	DiffEqRHS& DiffEqRHS::operator=(const DiffEqRHS& other) {
		if (this != &other)
		{
			_clean_up();

			Grid::operator=(other);
			_copy(other);
		};
		return *this;
	};
	DiffEqRHS& DiffEqRHS::operator=(DiffEqRHS&& other) {
		if (this != &other)
		{
			_clean_up();
			Grid::operator=(other);
			_move(other);
		};
		return *this;
	};

	DiffEqRHS::~DiffEqRHS() {
		_clean_up();
	};
	void DiffEqRHS::_copy(const DiffEqRHS& other) {
		_name = other._name;
		_no_dims = other._no_dims;
		_domain = other._domain;
		_parent_ixn_param = other._parent_ixn_param;
		_updates = other._updates;
	};
	void DiffEqRHS::_move(DiffEqRHS& other) {
		_name = other._name;
		_no_dims = other._no_dims;
		_domain = other._domain;
		_parent_ixn_param = other._parent_ixn_param;
		_updates = other._updates;

		// Reset other
		other._name = "";
		other._no_dims = 0;
		other._domain.clear();
		other._parent_ixn_param = nullptr;
		other._updates.clear();
	};

	void DiffEqRHS::_clean_up() {
		// Note: Domain is not owned => do not delete
	};

	/********************
	Validate
	********************/

	// Validate setup
	void DiffEqRHS::validate_setup() const {
		std::cout << "--- Validate Basis func " << _name << " ---" << std::endl;
	};

	/********************
	Getters
	********************/

	std::string DiffEqRHS::get_name() const {
		return _name;
	};

	Iptr DiffEqRHS::get_parent_ixn_param() const {
		return _parent_ixn_param;
	};

	const std::vector<Domain1D*>& DiffEqRHS::get_domain() const {
		return _domain;
	};

	std::vector<double> DiffEqRHS::_form_abscissas(int timepoint) const {
		std::vector<double> abscissas;
		for (auto dim=0; dim<_domain.size(); dim++) {
			abscissas.push_back(_domain[dim]->get_ixn_param()->get_val_at_timepoint(timepoint));
		};
		return abscissas;
	};

	bool DiffEqRHS::check_val_is_in_domain_at_timepoint(int timepoint) const {
		return Grid::check_in_domain(_form_abscissas(timepoint));
	};
	double DiffEqRHS::get_val_at_timepoint(int timepoint) const {
		return Grid::get_val(_form_abscissas(timepoint));
	};
	double DiffEqRHS::get_deriv_wrt_u_at_timepoint(int timepoint, q3c1::IdxSet global_vertex_idxs, std::vector<q3c1::DimType> dim_types) const {
		return Grid::get_deriv_wrt_coeff(_form_abscissas(timepoint), global_vertex_idxs, dim_types);
	};
	double DiffEqRHS::get_deriv_wrt_nu_at_timepoint(int timepoint, int deriv_dim) const {
		return Grid::get_deriv_wrt_abscissa(_form_abscissas(timepoint), deriv_dim);
	};

	/********************
	Update
	********************/

	void DiffEqRHS::update_calculate_and_store(int timepoint_start, int timepoint_end, double dt, double dopt) {

		// Adjoint
		double adjoint_val;

		// Abscissas
		std::vector<double> abscissas;

		// Updates
		std::map<q3c1::Vertex*,std::vector<double>> updates;

		// Go through all times
		for (auto timepoint=timepoint_start; timepoint<timepoint_end; timepoint++) {

			// Adjoint
			adjoint_val = _parent_ixn_param->get_adjoint()->get_val_at_timepoint(timepoint);

			// Form abscissas
			abscissas = _form_abscissas(timepoint);

			// Get updates for verts
			updates = Grid::get_deriv_wrt_coeffs_for_all_surrounding_verts(abscissas);

			// Append to any existing
			for (auto &pr: updates) {
				auto it = _updates.find(pr.first);
				if (it == _updates.end()) {
					// new
					_updates[pr.first] = pr.second;
					for (auto i=0; i<pr.second.size(); i++) {
						_updates[pr.first][i] *= dt * dopt * adjoint_val;
					};
				} else {
					// append
					for (auto i=0; i<pr.second.size(); i++) {
						_updates[pr.first][i] += dt * dopt * adjoint_val * pr.second[i];
					};
				};
			};
		};
	};

	// Verbose
	void DiffEqRHS::print_update_stored() const {
		for (auto const &pr: _updates) {
			std::cout << pr.first << ": " << std::flush;
			for (auto x: pr.second) {
				std::cout << x << " " << std::flush;
			};
			std::cout << std::endl;
		};
	};

	// Committ the update
	void DiffEqRHS::update_committ_stored(int opt_step, bool nesterov_mode) {

		// Go through stored updates
		if (nesterov_mode && opt_step > 1) {

			if (_no_dims == 1) {
				for (auto const &pr: _updates) {
					// Val
					pr.first->get_bf({q3c1::DimType::VAL})->increment_coeff((1.0 + (
						opt_step - 1.0) / (opt_step + 2.0)) * pr.second[0]);
					// Deriv
					pr.first->get_bf({q3c1::DimType::DERIV})->increment_coeff((1.0 + (
						opt_step - 1.0) / (opt_step + 2.0)) * pr.second[1]);
				};
			} else if (_no_dims == 2) {
				for (auto const &pr: _updates) {
					// Val-val
					pr.first->get_bf({q3c1::DimType::VAL,q3c1::DimType::VAL})->increment_coeff((1.0 + (
						opt_step - 1.0) / (opt_step + 2.0)) * pr.second[0]);
					// Val-Deriv
					pr.first->get_bf({q3c1::DimType::VAL,q3c1::DimType::DERIV})->increment_coeff((1.0 + (
						opt_step - 1.0) / (opt_step + 2.0)) * pr.second[1]);
					// Deriv-val
					pr.first->get_bf({q3c1::DimType::DERIV,q3c1::DimType::VAL})->increment_coeff((1.0 + (
						opt_step - 1.0) / (opt_step + 2.0)) * pr.second[2]);
					// Deriv-deriv
					pr.first->get_bf({q3c1::DimType::DERIV,q3c1::DimType::DERIV})->increment_coeff((1.0 + (
						opt_step - 1.0) / (opt_step + 2.0)) * pr.second[3]);
				};
			} else if (_no_dims == 3) {
				for (auto const &pr: _updates) {
					// Val-val-val
					pr.first->get_bf({q3c1::DimType::VAL,q3c1::DimType::VAL,q3c1::DimType::VAL})->increment_coeff((1.0 + (
						opt_step - 1.0) / (opt_step + 2.0)) * pr.second[0]);
					// Val-val-deriv
					pr.first->get_bf({q3c1::DimType::VAL,q3c1::DimType::VAL,q3c1::DimType::DERIV})->increment_coeff((1.0 + (
						opt_step - 1.0) / (opt_step + 2.0)) * pr.second[1]);
					// Val-deriv-val
					pr.first->get_bf({q3c1::DimType::VAL,q3c1::DimType::DERIV,q3c1::DimType::VAL})->increment_coeff((1.0 + (
						opt_step - 1.0) / (opt_step + 2.0)) * pr.second[2]);
					// Deriv-val-val
					pr.first->get_bf({q3c1::DimType::DERIV,q3c1::DimType::VAL,q3c1::DimType::VAL})->increment_coeff((1.0 + (
						opt_step - 1.0) / (opt_step + 2.0)) * pr.second[3]);
					// Val-deriv-deriv
					pr.first->get_bf({q3c1::DimType::VAL,q3c1::DimType::DERIV,q3c1::DimType::DERIV})->increment_coeff((1.0 + (
						opt_step - 1.0) / (opt_step + 2.0)) * pr.second[4]);
					// Deriv-val-deriv
					pr.first->get_bf({q3c1::DimType::DERIV,q3c1::DimType::VAL,q3c1::DimType::DERIV})->increment_coeff((1.0 + (
						opt_step - 1.0) / (opt_step + 2.0)) * pr.second[5]);
					// Deriv-deriv-val
					pr.first->get_bf({q3c1::DimType::DERIV,q3c1::DimType::DERIV,q3c1::DimType::VAL})->increment_coeff((1.0 + (
						opt_step - 1.0) / (opt_step + 2.0)) * pr.second[6]);
					// Deriv-deriv-deriv
					pr.first->get_bf({q3c1::DimType::DERIV,q3c1::DimType::DERIV,q3c1::DimType::DERIV})->increment_coeff((1.0 + (
						opt_step - 1.0) / (opt_step + 2.0)) * pr.second[7]);
				};
			};

		} else {

			if (_no_dims == 1) {
				for (auto const &pr: _updates) {
					// Val
					pr.first->get_bf({q3c1::DimType::VAL})->increment_coeff(pr.second[0]);
					// Deriv
					pr.first->get_bf({q3c1::DimType::DERIV})->increment_coeff(pr.second[1]);
				};
			} else if (_no_dims == 2) {
				for (auto const &pr: _updates) {
					// Val-val
					pr.first->get_bf({q3c1::DimType::VAL,q3c1::DimType::VAL})->increment_coeff(pr.second[0]);
					// Val-Deriv
					pr.first->get_bf({q3c1::DimType::VAL,q3c1::DimType::DERIV})->increment_coeff(pr.second[1]);
					// Deriv-val
					pr.first->get_bf({q3c1::DimType::DERIV,q3c1::DimType::VAL})->increment_coeff(pr.second[2]);
					// Deriv-deriv
					pr.first->get_bf({q3c1::DimType::DERIV,q3c1::DimType::DERIV})->increment_coeff(pr.second[3]);
				};
			} else if (_no_dims == 3) {
				for (auto const &pr: _updates) {
					// Val-val-val
					pr.first->get_bf({q3c1::DimType::VAL,q3c1::DimType::VAL,q3c1::DimType::VAL})->increment_coeff(pr.second[0]);
					// Val-val-deriv
					pr.first->get_bf({q3c1::DimType::VAL,q3c1::DimType::VAL,q3c1::DimType::DERIV})->increment_coeff(pr.second[1]);
					// Val-deriv-val
					pr.first->get_bf({q3c1::DimType::VAL,q3c1::DimType::DERIV,q3c1::DimType::VAL})->increment_coeff(pr.second[2]);
					// Deriv-val-val
					pr.first->get_bf({q3c1::DimType::DERIV,q3c1::DimType::VAL,q3c1::DimType::VAL})->increment_coeff(pr.second[3]);
					// Val-deriv-deriv
					pr.first->get_bf({q3c1::DimType::VAL,q3c1::DimType::DERIV,q3c1::DimType::DERIV})->increment_coeff(pr.second[4]);
					// Deriv-val-deriv
					pr.first->get_bf({q3c1::DimType::DERIV,q3c1::DimType::VAL,q3c1::DimType::DERIV})->increment_coeff(pr.second[5]);
					// Deriv-deriv-val
					pr.first->get_bf({q3c1::DimType::DERIV,q3c1::DimType::DERIV,q3c1::DimType::VAL})->increment_coeff(pr.second[6]);
					// Deriv-deriv-deriv
					pr.first->get_bf({q3c1::DimType::DERIV,q3c1::DimType::DERIV,q3c1::DimType::DERIV})->increment_coeff(pr.second[7]);
				};
			};
		};

		// Clear keys
		_updates.clear();

	};
};