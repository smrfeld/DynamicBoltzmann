#include "var_term_traj.hpp"
#include <iostream>
#include "general.hpp"
#include <fstream>

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Variational Term Trajectory
	****************************************/

	/********************
	Constructor
	********************/

	VarTermTraj::VarTermTraj(std::string name, IxnParamTraj *num, BasisFunc *denom, std::vector<IxnParamTraj*> denom_ixn_params, BasisFunc *num_bf, int n_ixn_params_in_num_bf, int n_t)
	{
		// Name
		_name = name;

		// Length (time)
		_n_t = n_t;

		// Num/denom ptrs
		_num = num;
		_denom = denom;

		// Check if it is a delta source
		_delta_source = _num->is_bf(denom);

		// Basis func corresponding to numerator
		_num_bf = num_bf;

		// Number interaction parameters in numerators basis func
		_n_ixn_params_in_num_bf=n_ixn_params_in_num_bf;

		// Derivatives of the num bf
		_num_bf_derivs = new double[_n_ixn_params_in_num_bf];

		// Vals
		for (int it=0; it<_n_t; it++) {
			_vals.push_back(Array(denom_ixn_params));
		};

		// Length of each array
		_val_len = 1;
		for (auto v: denom_ixn_params) { _val_len *= v->n(); };
	};
	VarTermTraj::VarTermTraj(const VarTermTraj& other)
	{
		_copy(other);
	};
	VarTermTraj& VarTermTraj::operator=(const VarTermTraj& other)
	{
		if (this != &other)
		{
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	VarTermTraj::~VarTermTraj()
	{
		_clean_up();
	};
	void VarTermTraj::_copy(const VarTermTraj& other) {
		_name = other._name;
		_n_t = other._n_t;
		_num = other._num;
		_denom = other._denom;
		_delta_source = other._delta_source;
		_num_bf = other._num_bf;
		_n_ixn_params_in_num_bf = other._n_ixn_params_in_num_bf;

		_num_bf_derivs = new double[_n_ixn_params_in_num_bf];
		std::copy(other._num_bf_derivs,other._num_bf_derivs+_n_ixn_params_in_num_bf,_num_bf_derivs);
		
		_vals = other._vals;
		_val_len = other._val_len;
		_update_var_terms = other._update_var_terms;
	};
	void VarTermTraj::_clean_up() {
		safeDelArr(_num_bf_derivs);
	};

	/********************
	Set the pointers needed to update this term
	********************/

	void VarTermTraj::add_update_ptr(VarTermTraj* var_term)
	{
		_update_var_terms.push_back(var_term);
	};

	/********************
	Validation
	********************/

	void VarTermTraj::validate_setup() const {
		std::cout << "--- Var term traj: " << _name << " ---" << std::endl;
		if (_num_bf) {
			std::cout << "   Numerator's basis func: " << _num_bf->name() << std::endl;
		} else {
			std::cerr << "ERROR: Var term traj: " << _name << " has no numerator basis func" << std::endl;
			exit(EXIT_FAILURE);
		};
		if (_delta_source) {
			std::cout << "   This is a delta source" << std::endl;
		};
		if (_update_var_terms.size() == 0) {
			std::cerr << "ERROR: No variational terms for updating" << std::endl;
			exit(EXIT_FAILURE);
		};
		for (auto vt_ptr: _update_var_terms) {
			std::cout << "   Updated using var term: " << vt_ptr->name() << std::endl;
		};
	};

	/********************
	// Calculate next timestep
	********************/

	void VarTermTraj::calculate_at_time(int it_next, double dt)
	{
		// Calculate derivative of basis funcs
		for (int j=0; j<_n_ixn_params_in_num_bf; j++) {
			_num_bf_derivs[j] = _num_bf->get_deriv_at_time(it_next-1, j);
		};

		// New val
		double d;

		// Iterate over all pts
		if (_delta_source)
		{
			for (int i=0; i<_val_len; i++) {
				d = _num_bf->get_delta_source(it_next-1, i);
				for (int j=0; j<_n_ixn_params_in_num_bf; j++) {
					d += _num_bf_derivs[j] * _update_var_terms[j]->get_at_time_by_idx(it_next-1, i);
				};
				// Set the value
				_vals[it_next].set_by_idx(i,_vals[it_next-1].get_by_idx(i)+dt*d);
			};
		} else {
			for (int i=0; i<_val_len; i++) {
				d = 0.0;
				for (int j=0; j<_n_ixn_params_in_num_bf; j++) {
					d += _num_bf_derivs[j] * _update_var_terms[j]->get_at_time_by_idx(it_next-1, i);
				};
				// Set the value
				_vals[it_next].set_by_idx(i,_vals[it_next-1].get_by_idx(i)+dt*d);
			};
		};
	};

	/********************
	Getters/setters
	********************/

	double VarTermTraj::get_at_time_by_idx(int it, int i) {
		return _vals[it].get_by_idx(i);
	};

	std::string VarTermTraj::name() {
		return _name;
	};

	/********************
	Write
	********************/

	void VarTermTraj::write_vals(std::string dir,int idx) const {
		std::ofstream f;
		f.open(dir+_name+"_"+pad_str(idx,4)+".txt");
		for (int i=0; i<_val_len; i++) {
			for (int t=0; t<_n_t; t++) {
				f << _vals[t].get_by_idx(i);
				if (t != _n_t-1) { f << " "; };
			};
			if (i!=_val_len-1) { f << "\n"; };
		};
		f.close();
	}; 

};