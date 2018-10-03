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

	DiffEqRHS::DiffEqRHS(std::string name, Iptr parent_ixn_param, std::vector<Domain1D*> domain) : dcu::Grid(std::vector<dcu::Dimension1D*>(domain.begin(),domain.end())), _idxs_k(domain.size()) {
		_name = name;
		_domain = domain;
		_no_dims = _domain.size();
		_parent_ixn_param = parent_ixn_param;

		// Init structures for evaluating
		_abscissas = new double[_no_dims];
		std::fill_n(_abscissas,_no_dims,0.0);
	};
	DiffEqRHS::DiffEqRHS(const DiffEqRHS& other) : Grid(other), _idxs_k(other._idxs_k) {
		_copy(other);
	};
	DiffEqRHS::DiffEqRHS(DiffEqRHS&& other) : Grid(std::move(other)), _idxs_k(std::move(other._idxs_k)) {
		_move(other);
	};

	DiffEqRHS& DiffEqRHS::operator=(const DiffEqRHS& other) {
		if (this != &other)
		{
			_clean_up();

			Grid::operator=(other);
			_idxs_k = other._idxs_k;
			_copy(other);
		};
		return *this;
	};
	DiffEqRHS& DiffEqRHS::operator=(DiffEqRHS&& other) {
		if (this != &other)
		{
			_clean_up();

			Grid::operator=(other);
			_idxs_k = std::move(other._idxs_k);
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
		_abscissas = new double[_no_dims];
		std::copy(other._abscissas,other._abscissas+_no_dims,_abscissas);
		_idxs_k = other._idxs_k;
		_parent_ixn_param = other._parent_ixn_param;
		_updates = other._updates;
	};
	void DiffEqRHS::_move(DiffEqRHS& other) {
		_name = other._name;
		_no_dims = other._no_dims;
		_domain = other._domain;
		_abscissas = other._abscissas;
		_idxs_k = other._idxs_k;
		_parent_ixn_param = other._parent_ixn_param;
		_updates = other._updates;

		// Reset other
		other._name = "";
		other._no_dims = 0;
		other._domain.clear();
		other._abscissas = nullptr;
		other._parent_ixn_param = nullptr;
		other._updates.clear();
	};

	void DiffEqRHS::_clean_up() {
		// Domain is not owned

		if (_abscissas) {
			delete[] _abscissas;
			_abscissas = nullptr;
		};
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

	void DiffEqRHS::_form_abscissas(int timepoint) {
		for (auto dim=0; dim<_domain.size(); dim++) {
			_abscissas[dim] = _domain[dim]->get_ixn_param()->get_val_at_timepoint(timepoint);
		};
	};

	double DiffEqRHS::get_val_at_timepoint(int timepoint) {
		_form_abscissas(timepoint);
		return get_val(_abscissas);
	};
	double DiffEqRHS::get_deriv_wrt_u_at_timepoint(int timepoint, dcu::IdxSet idxs_k) {
		_form_abscissas(timepoint);
		return get_deriv_wrt_pt_value(_abscissas, idxs_k);
	};
	double DiffEqRHS::get_deriv_wrt_nu_at_timepoint(int timepoint, int k) {
		_form_abscissas(timepoint);
		return get_deriv_wrt_x(_abscissas, k);
	};

	/********************
	Update
	********************/

	void DiffEqRHS::update_calculate_and_store(int timepoint_start, int timepoint_end, double dt, double dopt) {

		// Adjoint
		double adjoint_val;

		// Init idxs set
		dcu::IdxSet tmp(_no_dims);
		dcu::Nbr4 nbr4(tmp);

		// Derivative of F
		double f_deriv;

		// The update
		double update;

		// Gridptin
		dcu::GridPtIn* grid_pt=nullptr;

		// Iterator
		std::map<dcu::GridPtIn*,double>::iterator it;

		// Go through all times
		for (auto timepoint=timepoint_start; timepoint<timepoint_end; timepoint++) {

			// Adjoint
			adjoint_val = _parent_ixn_param->get_adjoint()->get_val_at_timepoint(timepoint);

			// Form abscissas
			_form_abscissas(timepoint);

			// Get surrounding nbrs
			nbr4 = get_surrounding_4_grid_pts(_abscissas);

			// Go through nbrs
			for (auto i=0; i<nbr4.get_no_grid_pts(); i++) {

				// Only update for grid pts that are inside!
				if (nbr4.get_grid_point(i)->get_type() == dcu::GridPtType::INSIDE) {

					// Get deriv
					f_deriv = get_deriv_wrt_pt_value(nbr4, nbr4.convert_idx_set(i));

					// The update
					update = dt * dopt * f_deriv * adjoint_val;

					// Grid pt
					grid_pt = nbr4.get_grid_point_inside(i);

					// Store
					it = _updates.find(grid_pt);
					if (it != _updates.end()) {
						_updates[grid_pt] -= update;
					} else {
						_updates[grid_pt] = -1.0 * update;
					};
				};
			};
		};
	};

	// Verbose
	void DiffEqRHS::print_update_stored() const {
		for (auto const &pr: _updates) {
			std::cout << *pr.first << ": " << pr.second << std::endl;
		};
	};

	// Committ the update
	void DiffEqRHS::update_committ_stored(int opt_step, bool nesterov_mode) {

		// Go through stored updates
		for (auto const &pr: _updates) {

			// If nesterov, move a bit further
			if (nesterov_mode && opt_step > 1) { // only start at opt step 2
				pr.first->increment_ordinate((1.0 + (
					opt_step - 1.0) / (opt_step + 2.0)) * pr.second);
			} else {
				pr.first->increment_ordinate(pr.second);
			};
		};

		// Clear keys
		_updates.clear();

	};
};