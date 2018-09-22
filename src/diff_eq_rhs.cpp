#include "../include/dynamicboltz_bits/diff_eq_rhs.hpp"

// Other headers
#include "../include/dynamicboltz_bits/general.hpp"
#include "../include/dynamicboltz_bits/ixn_param.hpp"

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

	DiffEqRHS::DiffEqRHS(std::string name, Domain domain) : dcu::Grid(domain.get_dimensions()), _idxs_k(domain.size()) {
		_name = name;
		_domain = domain.get_domain();

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
		_abscissas = other._abscissas;
	};
	void DiffEqRHS::_move(DiffEqRHS& other) {
		_name = other._name;
		_domain = other._domain;
		_abscissas = other._abscissas;
	};

	void DiffEqRHS::_clean_up() {};

	/********************
	Move to the nesterov intermediate point
	********************/
	/*
	void DiffEqRHS::nesterov_move_to_intermediate_pt(int opt_step) {
		if (!_nesterov_prev_pt) {
			std::cerr << "Error! No prev nesterov pt exists in basis func " << _name << std::endl;
			exit(EXIT_FAILURE);
		};

		// Move to the intermediate point
		double curr;
		for (int i=0; i<_val_len; i++) {
			// Tem store the old
			curr = _vals[i];

			// Move
			_vals[i] = curr + (opt_step - 1.0) / (opt_step + 2.0) * (curr - _nesterov_prev_pt->get_by_idx(i));

			// The current point will next be the old
			_nesterov_prev_pt->set_by_idx(i, curr);
		};
	};
	*/
	
	/********************
	Set prev nesterov
	********************/

	/*
	void DiffEqRHS::nesterov_set_prev_equal_curr() {
		if (!_nesterov_prev_pt) {
			// Make
			_nesterov_prev_pt = new Array(_domain);
		};
		// Copy
		for (int i=0; i<_val_len; i++) {
			_nesterov_prev_pt->set_by_idx(i, _vals[i]);
		};
	};
	*/

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
	Calculate the new basis function
	********************/

	/*
	void DiffEqRHS::update(int t_start, int t_end, double dt, double dopt, bool exp_decay, double exp_decay_t0, double exp_decay_lambda, bool l2_reg_params_mode, std::map<IxnParam*,double> &l2_lambda_params, std::map<IxnParam*,double> &l2_reg_centers) 
	{
		int *idxs;
		double *nu_vals;
		double decay = 1.0;

		// Go through all idxs
		double up1, up2, l2_center;
		for (int i=0; i<_val_len; i++) {

			// Go through all times
			for (int t=t_start; t<t_end; t++) {

				// Exp decay
				if (exp_decay) {
					decay = exp(-exp_decay_lambda*abs(t - exp_decay_t0));
				};
	*/
				// Go through all updating terms
				/*
				for (auto p: _update_ptrs) {
					up1 = dopt * dt * p.first->moments_diff_at_time(t) * p.second->get_val_at_timepoint_by_idx(t, i) * decay;
					_vals[i] += up1;

					// L2
					if (l2_reg_params_mode) {
						// Get numerator of var term
						IxnParam* num = p.second->get_numerator_ixn_param_traj();

						// Lookup l2 lambda
						auto it = l2_lambda_params.find(num);
						if (it != l2_lambda_params.end()) {
							// Lookup center if it exists, else 0
							auto it2 = l2_reg_centers.find(num);
							if (it2 != l2_reg_centers.end()) {
								l2_center = it2->second;
							} else {
								l2_center = 0.;
							};
							up2 = dopt * it->second * dt * sgn(num->get_val_at_timepoint(t) - l2_center) * abs(num->get_val_at_timepoint(t) - l2_center) * p.second->get_val_at_timepoint_by_idx(t, i) * decay;

							_vals[i] -= up2;
							// std::cout << up1 << " " << up2 << std::endl;
						};
					};
				};
				*/
	/*
			};
		};
	};
	*/

	/*
	void DiffEqRHS::update_gather(int t_start, int t_end, double dt, double dopt, bool exp_decay, double exp_decay_t0, double exp_decay_lambda, bool l2_reg_params_mode, std::map<IxnParam*,double> &l2_lambda_params, std::map<IxnParam*,double> &l2_reg_centers) 
	{
		if (!_update_gathered) {
			// alloc
			_update_gathered = new double[_val_len];
			std::fill_n(_update_gathered,_val_len,0.);
		};

		int *idxs;
		double *nu_vals;
		double decay = 1.0;
		double up1,up2,l2_center;
		IxnParam* num;

		// Go through all idxs
		for (int i=0; i<_val_len; i++) {

			// Go through all times
			for (int t=t_start; t<t_end; t++) {

				// Exp decay
				if (exp_decay) {
					decay = exp(-exp_decay_lambda*abs(t - exp_decay_t0));
				};
	*/
				// Go through all updating terms
				/*
				for (auto p: _update_ptrs) {
					up1 = p.first->moments_diff_at_time(t);

					// L2
					if (l2_reg_params_mode) {
						// Get numerator of var term
						num = p.second->get_numerator_ixn_param_traj();

						// Lookup l2 lambda
						auto it = l2_lambda_params.find(num);
						if (it != l2_lambda_params.end()) {
							// Lookup center if it exists, else 0
							auto it2 = l2_reg_centers.find(num);
							if (it2 != l2_reg_centers.end()) {
								l2_center = it2->second;
							} else {
								l2_center = 0.;
							};
							up2 = it->second * sgn(num->get_val_at_timepoint(t) - l2_center) * abs(num->get_val_at_timepoint(t) - l2_center);

							// Subtract from up1
							up1 -= up2;
						};
					};

					_update_gathered[i] += dopt * dt * up1 * p.second->get_val_at_timepoint_by_idx(t, i) * decay / (t_end - t_start);
				};
				*/
	/*
			};
		};
	};

	void DiffEqRHS::update_committ_gathered() 
	{
		if (!_update_gathered) {
			std::cerr << "ERROR! No update allocated." << std::endl;
			exit(EXIT_FAILURE);
		};

		// Go through all idxs
		for (int i=0; i<_val_len; i++) {
			_vals[i] += _update_gathered[i];
		};

		// Reset to 0
		std::fill_n(_update_gathered,_val_len,0.);
	};
	*/


};