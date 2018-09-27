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
	Domain
	****************************************/

	Domain::Domain(std::vector<Domain1D> domain) {
		_domain = domain;
		for (auto &d: _domain) {
			_dimensions.push_back(dcu::Dimension1D(d));
		};
	};
	Domain::Domain(const Domain& other) {
		_copy(other);
	};
	Domain::Domain(Domain&& other) {
		_copy(other);
		other._reset();
	};
	Domain& Domain::operator=(const Domain& other) {
		if (this != &other)
		{
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	Domain& Domain::operator=(Domain&& other) {
		if (this != &other)
		{
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
	};
	Domain::~Domain() {
		_clean_up();
	};
	void Domain::_copy(const Domain& other)
	{
		_dimensions = other._dimensions;
		_domain = other._domain;
	};
	void Domain::_reset()
	{
		_dimensions.clear();
		_domain.clear();
	};
	void Domain::_clean_up() {
	};

	/********************
	Setters
	********************/

	void Domain::add_dimension(Domain1D domain) {
		_domain.push_back(domain);
		_dimensions.push_back(dcu::Dimension1D(domain));
	};

	/********************
	Getters
	********************/

	int Domain::size() const {
		return _domain.size();
	};
	std::vector<dcu::Dimension1D> Domain::get_dimensions() const {
		return _dimensions;
	};
	std::vector<Domain1D> Domain::get_domain() const {
		return _domain;
	};









































	/****************************************
	DiffEqRHS
	****************************************/

	/********************
	Constructor
	********************/

	DiffEqRHS::DiffEqRHS(std::string name, Iptr parent_ixn_param, Domain domain) : dcu::Grid(domain.get_dimensions()), _idxs_k(domain.size()) {
		_name = name;
		_domain = domain.get_domain();
		_dimensions = domain.get_dimensions();
		_parent_ixn_param = parent_ixn_param;

		// Init structures for evaluating
		for (auto dim=0; dim<_domain.size(); dim++) {
			_abscissas.push_back(0.0);
		};
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
		_domain = other._domain;
		_dimensions = other._dimensions;
		_abscissas = other._abscissas;
		_parent_ixn_param = other._parent_ixn_param;
		_updates = other._updates;
	};
	void DiffEqRHS::_move(DiffEqRHS& other) {
		_name = other._name;
		_domain = other._domain;
		_dimensions = other._dimensions;
		_abscissas = other._abscissas;
		_parent_ixn_param = other._parent_ixn_param;
		_updates = other._updates;
	};

	void DiffEqRHS::_clean_up() {};

	/********************
	Nesterov
	********************/

	// Move to the nesterov intermediate point
	void DiffEqRHS::nesterov_move_to_intermediate_pt(int i_opt_step) {
		// Only start at step 2; else intermediate pt is current
		if (i_opt_step <= 1) {
			return;
		};

		// Go through updates of prev step
		// current - prev = current - (current - update) = update
		for (auto &pr: _updates) {
			get_grid_point_ref(pr.first).increment_ordinate(
				(i_opt_step - 1.0) / (i_opt_step + 2.0) * pr.second
				);
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

	const std::vector<Domain1D>& DiffEqRHS::get_domain() const {
		return _domain;
	};

	void DiffEqRHS::_form_abscissas(int timepoint) {
		for (auto dim=0; dim<_domain.size(); dim++) {
			_abscissas[dim] = _domain[dim].get_ixn_param()->get_val_at_timepoint(timepoint);
		};
	};

	double DiffEqRHS::get_val_at_timepoint(int timepoint) {
		_form_abscissas(timepoint);
		return get_val_by_ref(_abscissas);
	};
	double DiffEqRHS::get_deriv_wrt_u_at_timepoint(int timepoint, dcu::IdxSet4 idxs_k) {
		_form_abscissas(timepoint);
		return get_deriv_wrt_pt_value_by_ref(_abscissas, idxs_k);
	};
	double DiffEqRHS::get_deriv_wrt_nu_at_timepoint(int timepoint, int k) {
		_form_abscissas(timepoint);
		return get_deriv_wrt_x_by_ref(_abscissas, k);
	};

	/********************
	Update
	********************/

	void DiffEqRHS::update_calculate_and_store(int timepoint_start, int timepoint_end, double dt, double dopt, bool clear_existing_updates) {

		// Clear previous updates
		if (clear_existing_updates) {
			_updates.clear();
		};

		// Adjoint
		double adjoint_val;

		// Init idxs set
		dcu::IdxSet idxs_global(_domain.size());
		dcu::IdxSet4 idxs_4(_domain.size());
		dcu::Nbr4 nbr4(idxs_global);

		// Derivative of F
		double f_deriv;

		// Go through all times
		for (auto timepoint=timepoint_start; timepoint<timepoint_end; timepoint++) {

			// Adjoint
			adjoint_val = _parent_ixn_param->get_adjoint()->get_val_at_timepoint(timepoint);

			// Form abscissas
			_form_abscissas(timepoint);

			// Get surrounding nbrs
			nbr4 = get_surrounding_4_grid_pts_by_ref(_abscissas);

			// Go through nbrs
			for (auto const &pr: nbr4.types) {

				// Only update for grid pts that are inside!
				if (pr.second == dcu::GridPtType::INSIDE) {

					// Get local idxs
					idxs_4 = pr.first;

					// Get deriv
					f_deriv = get_deriv_wrt_pt_value_by_ref(_abscissas, idxs_4);

					// Convert to global idxs
					// 0 = i-1
					idxs_global = nbr4.idxs_i - 1 + idxs_4;

					// Store update
					_updates[dcu::GridPtKey(idxs_global,_dimensions)] -= dt * dopt * f_deriv * adjoint_val;

				};
			};
		};
	};

	// Verbose
	void DiffEqRHS::print_update_stored() const {
		for (auto const &pr: _updates) {
			std::cout << "Stored update for grid pt: " << pr.first << " = " << pr.second << std::endl;
		};
	};

	// Committ the update
	void DiffEqRHS::update_committ_stored() {

		// Go through stored updates
		for (auto const &pr: _updates) {

			// Increment
			get_grid_point_ref(pr.first).increment_ordinate(pr.second);

		};

	};
};