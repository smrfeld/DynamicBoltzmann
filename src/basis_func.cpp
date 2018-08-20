#include "basis_func.hpp"

// Other headers
#include "../include/dynamicboltz_bits/general.hpp"
#include "ixn_param_traj.hpp"
#include "var_term_traj.hpp"

#include <iostream>
#include "math.h"
#include <fstream>
#include <iomanip>

/************************************
* Namespace for dboltz
************************************/

namespace dboltz {

	// Sign function
	int sgn(double val) {
		return (0. < val) - (val < 0.);
	};

	/****************************************
	Array
	****************************************/

	Array::Array(std::vector<IxnParamTraj*> ixn_params)
	{
		_ixn_params = ixn_params;
		_n_params = _ixn_params.size();

		// Values
		_val_len = 1;
		for (auto v: _ixn_params) { _val_len *= v->n(); };
		_vals = new double[_val_len];

		// Zero by default
		std::fill_n(_vals,_val_len,0.0);

		// Update dimension powers
		int pwr;
		for (int iv=0; iv<_n_params; iv++) 
		{
			pwr = 1;
			for (int jv=iv+1; jv<_n_params; jv++)
			{
				pwr *= _ixn_params[jv]->n();
			};
			_dim_pwrs.push_back(pwr); 
		};
	};
	Array::Array(IxnParamTraj* ixn_param) : Array(std::vector<IxnParamTraj*>{ixn_param}) {};
	Array::Array(const Array& other) 
	{
		_copy(other);
	};
	Array& Array::operator=(const Array& other)
	{
		if (this != &other) {
			_clean_up();
			_copy(other);	
		};
		return *this;
	};
	Array::Array(Array&& other)
	{
		// steal other's resources
		_n_params = other._n_params;
		_val_len = other._val_len;
		_ixn_params = other._ixn_params;
		_dim_pwrs = other._dim_pwrs;
		_vals = new double[_val_len];
		std::copy( other._vals, other._vals + _val_len, _vals );
		// reset other
		other._n_params = 0;
		other._val_len = 0;
		other._ixn_params.clear();
		other._dim_pwrs.clear();
		safeDelArr(other._vals);
	};
	Array& Array::operator=(Array&& other)
	{
		if (this!=&other)
		{
			// release the current object’s resources
			_clean_up();
			// steal other’s resource
			_n_params = other._n_params;
			_val_len = other._val_len;
			_ixn_params = other._ixn_params;
			_dim_pwrs = other._dim_pwrs;
			_vals = new double[_val_len];
			std::copy( other._vals, other._vals + _val_len, _vals );
			// reset other
			other._n_params = 0;
			other._val_len = 0;
			other._ixn_params.clear();
			other._dim_pwrs.clear();
			safeDelArr(other._vals);		
		};
		return *this;
	};
	Array::~Array()
	{
		_clean_up();
	};
	void Array::_copy(const Array& other) {
		_n_params = other._n_params;
		_val_len = other._val_len;
		_ixn_params = other._ixn_params;
		_dim_pwrs = other._dim_pwrs;

		_vals = new double[_val_len];
		std::copy( other._vals, other._vals + _val_len, _vals );
	};
	void Array::_clean_up() {
		safeDelArr(_vals);
	};

	/********************
	Get/set an element by index
	********************/

	double Array::get_by_idxs(int *idxs) const
	{
		int i=0;
		for (int id=0; id<_n_params; id++) {
			i += _dim_pwrs[id] * idxs[id];
		};
		return _vals[i];
	};

	double Array::get_by_idx(int i) const
	{
		return _vals[i];
	};

	void Array::set_by_idxs(int *idxs, double val)
	{
		int i=0;
		for (int id=0; id<_n_params; id++) {
			i += _dim_pwrs[id] * idxs[id];
		};
		_vals[i] = val;
	};

	void Array::set_by_idx(int i, double val)
	{
		_vals[i] = val;
	};

	/********************
	Get indexes by element
	********************/

	void Array::get_idxs(int i, int *idxs) const {
		int rem = i;
		for (int id=0; id<_n_params; id++) {
			idxs[id] = int(rem/_dim_pwrs[id]);
			rem = rem % _dim_pwrs[id];
		};
	};

	/********************
	Write to file
	********************/

	void Array::write_grid(std::string fname) const
	{
		std::ofstream f;
		f.open (fname);

		// Run through all points
		int *idxs = new int[_n_params];
		for (int i=0; i<_val_len; i++)
		{
			get_idxs(i, idxs);
			for (int ip=0; ip<_n_params; ip++) 
			{
				f << _ixn_params[ip]->get_by_idx(idxs[ip]);
				if (ip != _n_params-1) { f << " "; };
			};
			f << "\n";
		};

		f.close();

		// Clean up
		safeDelArr(idxs);
	};

	void Array::write_vals(std::string dir, std::string name, int idx_opt_step) const
	{
		std::ofstream f;
		f.open (dir+name+"_"+pad_str(idx_opt_step,4)+".txt");
		for (int i=0; i<_val_len; i++)
		{
			f << std::setprecision(15) << _vals[i];
			if (i!=_val_len-1) { f << "\n"; };
		};
		f.close();
	};

	void Array::read_vals(std::string fname) 
	{
		std::ifstream f;
		f.open(fname);
		char frag[100]; // fragments of the line
		std::string x="";
		int i=0;
		if (f.is_open()) { // make sure we found it
			while (!f.eof()) {
			    f >> frag;
			    _vals[i] = atof(frag);
			    i++;
			};
		};
		f.close();
	};

	/********************
	Check dimensions against another array
	********************/

	bool Array::check_dims(const Array& other) const
	{
		// Check no dims
		if (_n_params != other._n_params) { 
			std::cerr << "ERROR! No dimensions don't match in update step" << std::endl;
			return false;
		};

		// Check each dim
		for (int id=0; id<_n_params; id++) {
			if (_ixn_params[id] != other._ixn_params[id]) {
				return false;
			};
		};
		return true;
	};

	/********************
	Zero
	********************/

	void Array::zero()
	{
		std::fill_n(_vals,_val_len,0.);
	};

	/****************************************
	BasisFunc
	****************************************/

	/********************
	Constructor
	********************/

	BasisFunc::BasisFunc(std::string name, IxnParamTraj* parent_ixn_param, std::vector<IxnParamTraj*> ixn_params) : Array(ixn_params) {
		_name = name;
		
		_parent_ixn_param = parent_ixn_param;

		_derivs = new bool[_n_params];
		std::fill_n(_derivs,_n_params,false);

		_idxs_bounding = new int[_n_params];
		_idxs_p_cube = new int[_n_params];
		_idxs_ext_1 = new int[_n_params];
		_idxs_ext_2 = new int[_n_params];
		_fracs = new double[_n_params];
		_p_cube = new double[int(pow(4,_n_params))];

		_update_gathered = nullptr; // allocated later if needed

		_nesterov_prev_pt = nullptr;
	};
	BasisFunc::BasisFunc(const BasisFunc& other) : Array(other) {
		_copy(other);
	};
	BasisFunc& BasisFunc::operator=(const BasisFunc& other) {
		if (this != &other)
		{
			Array::operator=(other);

			_clean_up();

			_copy(other);
		};
		return *this;
	};
	BasisFunc::~BasisFunc() {
		_clean_up();
	};
	void BasisFunc::_copy(const BasisFunc& other) {
		_name = other._name;
		_parent_ixn_param = other._parent_ixn_param;
		_update_ptrs = other._update_ptrs;
		_derivs = new bool[_n_params];
		_idxs_bounding = new int[_n_params];
		_idxs_p_cube = new int[_n_params];
		_idxs_ext_1 = new int[_n_params];
		_idxs_ext_2 = new int[_n_params];
		_fracs = new double[_n_params];
		_p_cube = new double[int(pow(4,_n_params))];

		std::copy( other._derivs, other._derivs + _n_params, _derivs );
		std::copy( other._idxs_bounding, other._idxs_bounding + _n_params, _idxs_bounding );
		std::copy( other._idxs_p_cube, other._idxs_p_cube + _n_params, _idxs_p_cube );
		std::copy( other._idxs_ext_1, other._idxs_ext_1 + _n_params, _idxs_ext_1 );
		std::copy( other._idxs_ext_2, other._idxs_ext_2 + _n_params, _idxs_ext_2 );
		std::copy( other._fracs, other._fracs + _n_params, _fracs );
		std::copy( other._p_cube, other._p_cube + int(pow(4,_n_params)), _p_cube );

		if (other._update_gathered) {
			_update_gathered = new double[_val_len];
			std::copy( other._update_gathered, other._update_gathered + _val_len, _update_gathered );
		} else {
			_update_gathered = nullptr;
		};

		if (other._nesterov_prev_pt) {
			_nesterov_prev_pt = new Array(*(other._nesterov_prev_pt));
		} else {
			_nesterov_prev_pt = nullptr;
		};
	};
	void BasisFunc::_clean_up() {
		safeDelArr(_derivs);
		safeDelArr(_idxs_bounding);
		safeDelArr(_idxs_p_cube);
		safeDelArr(_idxs_ext_1);
		safeDelArr(_idxs_ext_2);
		safeDelArr(_fracs);
		safeDelArr(_p_cube);
		safeDelArr(_update_gathered);
		if (_nesterov_prev_pt) {
			delete _nesterov_prev_pt;
			_nesterov_prev_pt = nullptr;
		};
	};

	/********************
	Move to the nesterov intermediate point
	********************/

	void BasisFunc::nesterov_move_to_intermediate_pt(int opt_step) {
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

	/********************
	Set prev nesterov
	********************/

	void BasisFunc::nesterov_set_prev_equal_curr() {
		if (!_nesterov_prev_pt) {
			// Make
			_nesterov_prev_pt = new Array(_ixn_params);
		};
		// Copy
		for (int i=0; i<_val_len; i++) {
			_nesterov_prev_pt->set_by_idx(i, _vals[i]);
		};
	};

	/********************
	Add pointers needed to update
	********************/

	void BasisFunc::add_update_ptrs(IxnParamTraj* ixn_param, VarTermTraj* var_term)
	{
		_update_ptrs.push_back(std::make_pair(ixn_param,var_term));
	};

	/********************
	Validate
	********************/

	// Validate setup
	void BasisFunc::validate_setup() const {
		std::cout << "--- Validate Basis func " << _name << " ---" << std::endl;
		if (_update_ptrs.size() == 0) {
			std::cerr << "ERROR: no update ptrs" << std::endl;
			exit(EXIT_FAILURE);
		};
		for (auto pr: _update_ptrs) {
			std::cout << "   Updated using ixn param: " << pr.first->name() << " and var term: " << pr.second->name() << std::endl;
		};
	};

	/********************
	Get values, if they are in the lattice
	********************/

	// Recursion to fill p
	// idx_deltas will in the inner most loop be -1, 0, 1, 2 for i-1,i,i+1,i+2
	// Initial dim parameter = 0
	// idxs = contains the i
	void BasisFunc::_fill_p(int dim)
	{
		// i-1 to i+2
		for(_idxs_p_cube[dim] = 0; _idxs_p_cube[dim] < 4; ++_idxs_p_cube[dim]) 
		{
			if (dim != _n_params-1)
			{
				// Inception - need another for loop
				_fill_p(dim+1);
			} else {
				// Do something!

				// Print
				/*
				std::cout << "p cube idxs: ";
				for (int a=0; a<_n_params; a++) {
					std::cout << _idxs_p_cube[a] << " ";
				};
				*/

				// Check if point is outside box
				bool outside=false; // Is it outside on at least one dim?
				int idx;
				for (int d=0; d<_n_params; d++) {
					idx = _idxs_bounding[d]+_idxs_p_cube[d]-1;
					if (idx < 0) {
						outside = true; // Yes it's outside on at least one dim
						_idxs_ext_1[d] = 0;
						_idxs_ext_2[d] = 1;
					} else if (idx > _ixn_params[d]->n()-1) {
						outside = true; // Yes it's outside on at least one dim
						_idxs_ext_1[d] = _ixn_params[d]->n()-1;
						_idxs_ext_2[d] = _ixn_params[d]->n()-2;
					} else {
						_idxs_ext_1[d] = idx;
						_idxs_ext_2[d] = idx;
					};
				};

				// Index in _p_cube to write to
				int i_write=0;
				for (int d=0; d<_n_params; d++) {
					i_write += pow(4,_n_params-d-1) * _idxs_p_cube[d];
				};

				// Finally, write
				if (outside) {
					_p_cube[i_write] = 2*get_by_idxs(_idxs_ext_1) - get_by_idxs(_idxs_ext_2);
					// std::cout << "writing: 2*(" << _idxs_ext_1[0] << "," << _idxs_ext_1[1] << ") - (" << _idxs_ext_2[0] << "," << _idxs_ext_2[1] << ")";
				} else {
					_p_cube[i_write] = get_by_idxs(_idxs_ext_1);
					// std::cout << "writing: (" << _idxs_ext_1[0] << "," << _idxs_ext_1[1] << ")";
				};
				// std::cout << std::endl;
			};
		};
	};

	// Function to get bounding cube of 4 points
	void BasisFunc::_get_bounding(int it, bool safe)
	{
		// Check that x is in the box
		if (safe) {
			for (int id=0; id<_n_params; id++) {
				if (!(_ixn_params[id]->in_grid(_ixn_params[id]->get_at_time(it)))) {
					std::cerr << "ERROR - " << _ixn_params[id]->get_at_time(it) << " is outside the grid:" << std::endl;
					_ixn_params[id]->print_grid_range();
					exit(EXIT_FAILURE);
				};
			};
		};

		// Get the interval - i-1,i,POINT,i+1,i+2
		// indexes in a dimension of length n run from 0 to n-1
		// this should return i = 0 to n-2
		// left boundary: i=0 -> -1,0,1,2 where point i=-1 needs extrapolation later
		// right boundary: i=n-2 -> n-3,n-2,n-1,n where point i=n needs extrapolation later
		// ALSO: get the fraction this point is between two successive, i.e. between i and i+1
		double x;
		for (int id=0; id<_n_params; id++) {
			x = _ixn_params[id]->get_at_time(it);
			_idxs_bounding[id] = _ixn_params[id]->surrounding_idxs(x);
			_fracs[id] = _ixn_params[id]->frac_between(x,_idxs_bounding[id]);
		};

		// Get bounding box
		_fill_p(0);
	};

	/********************
	Get values, if they are in the lattice
	********************/

	double BasisFunc::get_at_time(int it) 
	{
		// Create the bounding box
		_get_bounding(it);

		return _n_cubic_interp(_n_params,_p_cube,_fracs,_derivs);
	};

	double BasisFunc::get_deriv_at_time(int it, int i_dim) {
		// Deriv
		_derivs[i_dim] = true;

		// Create the bounding box
		_get_bounding(it);

		// Eval
		double x = _n_cubic_interp(_n_params,_p_cube,_fracs,_derivs);

		// Reset
		_derivs[i_dim] = false;

		return x;
	};

	/********************
	Name
	********************/

	std::string BasisFunc::name() const {
		return _name;
	};

	/********************
	Calculate the new basis function
	********************/

	void BasisFunc::update(int t_start, int t_end, double dt, double dopt, bool exp_decay, double exp_decay_t0, double exp_decay_lambda, bool l2_reg_params_mode, std::map<IxnParamTraj*,double> &l2_lambda_params, std::map<IxnParamTraj*,double> &l2_reg_centers) 
	{
		int *idxs;
		double *nu_vals;
		double decay = 1.0;

		// L2 for params
		/*
		double l2_params=0.0;
		if (l2_reg_params_mode) {
			l2_params = 0.;
			// Go through all times
			for (int t=t_start; t<t_end; t++) {
				// Get the l2
				auto f = l2_lambda_params.find(_parent_ixn_param);
				if (f != l2_lambda_params.end()) {
					l2_params += dopt * f->second * dt * sgn(_parent_ixn_param->get_at_time(t)) * abs(_parent_ixn_param->get_at_time(t));
				};
			};
		};
		*/

		// Go through all idxs
		double up1, up2, l2_center;
		for (int i=0; i<_val_len; i++) {

			// Go through all times
			for (int t=t_start; t<t_end; t++) {

				// Exp decay
				if (exp_decay) {
					decay = exp(-exp_decay_lambda*abs(t - exp_decay_t0));
				};

				// Go through all updating terms
				for (auto p: _update_ptrs) {
					up1 = dopt * dt * p.first->moments_diff_at_time(t) * p.second->get_at_time_by_idx(t, i) * decay;
					_vals[i] += up1;

					// L2
					if (l2_reg_params_mode) {
						// Get numerator of var term
						IxnParamTraj* num = p.second->get_numerator_ixn_param_traj();

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
							up2 = dopt * it->second * dt * sgn(num->get_at_time(t) - l2_center) * abs(num->get_at_time(t) - l2_center) * p.second->get_at_time_by_idx(t, i) * decay;

							_vals[i] -= up2;
							// std::cout << up1 << " " << up2 << std::endl;
						};
					};
				};
			};

			// L2 for params
			/*
			if (l2_reg_params_mode) {
				_vals[i] -= l2_params;
			};
			*/
		};
	};

	void BasisFunc::update_gather(int t_start, int t_end, double dt, double dopt, bool exp_decay, double exp_decay_t0, double exp_decay_lambda) 
	{
		if (!_update_gathered) {
			// alloc
			_update_gathered = new double[_val_len];
			std::fill_n(_update_gathered,_val_len,0.);
		};

		int *idxs;
		double *nu_vals;
		double decay = 1.0;

		// Go through all idxs
		for (int i=0; i<_val_len; i++) {

			// Go through all times
			for (int t=t_start; t<t_end; t++) {

				// Exp decay
				if (exp_decay) {
					decay = exp(-exp_decay_lambda*abs(t - exp_decay_t0));
				};

				// Go through all updating terms
				for (auto p: _update_ptrs) {
					_update_gathered[i] += dopt * dt * p.first->moments_diff_at_time(t) * p.second->get_at_time_by_idx(t, i) * decay;
				};
			};
		};
	};

	void BasisFunc::update_committ_gathered(int t_start, int t_end, double dt, double dopt, bool l2_reg_params_mode, std::map<IxnParamTraj*,double> &l2_lambda_params, std::map<IxnParamTraj*,double> &l2_reg_centers) 
	{
		if (!_update_gathered) {
			std::cerr << "ERROR! No update allocated." << std::endl;
			exit(EXIT_FAILURE);
		};

		double l2_params=0.0;

		// L2 for params
		if (l2_reg_params_mode) {
			std::cout << "WARNING: There may be an error here for the L2!" << std::endl;

			l2_params = 0.;
			// Go through all times
			for (int t=t_start; t<t_end; t++) {
				// Get the l2
				auto f = l2_lambda_params.find(_parent_ixn_param);
				if (f != l2_lambda_params.end()) {
					l2_params += dopt * f->second * dt * sgn(_parent_ixn_param->get_at_time(t)) * abs(_parent_ixn_param->get_at_time(t));
				};
			};
		};

		// Go through all idxs
		for (int i=0; i<_val_len; i++) {
			_vals[i] += _update_gathered[i];

			// L2 for params
			if (l2_reg_params_mode) {
				_vals[i] -= l2_params;
			};
		};

		// Reset to 0
		std::fill_n(_update_gathered,_val_len,0.);
	};

	/********************
	Test fill in various dimensions
	********************/

	void BasisFunc::test_fill_2d() {
		if (_n_params != 2) { return; };
		std::vector<double> x0 = _ixn_params[0]->test_sin();
		std::vector<double> x1 = _ixn_params[1]->test_cos();

		int idxs[2];
		for (int i=0; i<x0.size(); i++) {
			for (int j=0; j<x1.size(); j++) {
				idxs[0] = i;
				idxs[1] = j;
				set_by_idxs(idxs, x0[i]*x1[j]);
			};
		};
	};
	void BasisFunc::test_fill_3d() {
		if (_n_params != 3) { return; };
		std::vector<double> x0 = _ixn_params[0]->test_sin();
		std::vector<double> x1 = _ixn_params[1]->test_cos();
		std::vector<double> x2 = _ixn_params[2]->test_cos();

		int idxs[3];
		for (int i=0; i<x0.size(); i++) {
			for (int j=0; j<x1.size(); j++) {
				for (int k=0; k<x2.size(); k++) {
					idxs[0] = i;
					idxs[1] = j;
					idxs[2] = k;
					set_by_idxs(idxs, x0[i]*x1[j]*x2[k]);
				};
			};
		};
	};

	/********************
	From parent
	********************/

	double BasisFunc::get_by_idxs(int *idxs) const {
		return Array::get_by_idxs(idxs);
	};

	double BasisFunc::get_by_idx(int i) const {
		return Array::get_by_idx(i);
	};

	void BasisFunc::set_by_idxs(int *idxs, double val) {
		Array::set_by_idxs(idxs,val);
	};

	void BasisFunc::write_grid(std::string fname) const {
		Array::write_grid(fname);
	}; 
	void BasisFunc::write_vals(std::string dir,int idx) const {
		Array::write_vals(dir,_name,idx);
	}; 
	void BasisFunc::read_vals(std::string fname) {
		Array::read_vals(fname);
	};

	/********************
	Get delta source
	********************/

	double BasisFunc::get_delta_source(int it, int i)
	{
		// Get idxs for this i
		int *idxs = new int[_n_params];
		get_idxs(i, idxs);

		// Go through ixn parmas
		double r=1;
		for (int ip=0; ip<_n_params; ip++) {
			r *= exp(-pow(_ixn_params[ip]->get_at_time(it)-_ixn_params[ip]->get_by_idx(idxs[ip]),2)/(2.*1.0*_ixn_params[ip]->delta())) / sqrt(2.*M_PI*1.0*_ixn_params[ip]->delta());
		};

		// Clean
		safeDelArr(idxs);

		return r;
	};

	/****************************************
	BasisFunc - PRIVATE
	****************************************/

	/********************
	Interpolation
	********************/

	double BasisFunc::_cubic_interp(double p[4], double f) {
		return p[1] + 0.5 * f *(p[2] - p[0] + f*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + f*(3.0*(p[1] - p[2]) + p[3] - p[0])));
	};
	double BasisFunc::_deriv(double p[4], double f) {
		return 0.5 * (p[2] - p[0]) + 2.0 * f * (p[0] - 2.5 * p[1] + 2.0 * p[2] - 0.5 * p[3]) + 3.0 * f * f *( - 0.5 * p[0] + 1.5 * p[1] - 1.5 * p[2] + 0.5 * p[3]);
	};
	double BasisFunc::_n_cubic_interp(int dim, double* p, double fracs[], bool derivs[])
	{
		if (dim == 1)
		{
			// Finall reached correct dim for linear interp
			if (*derivs) {
				return _deriv(p, *fracs);
			} else {
 				return _cubic_interp(p, *fracs);
			};
		} else {
			// Interpolate in next dim
			double arr[4];
			int skip = 1 << (dim - 1) * 2;
			arr[0] = _n_cubic_interp(dim-1, p, fracs+1, derivs+1);
			arr[1] = _n_cubic_interp(dim-1, p + skip, fracs+1, derivs+1);
			arr[2] = _n_cubic_interp(dim-1, p + 2*skip, fracs+1, derivs+1);
			arr[3] = _n_cubic_interp(dim-1, p + 3*skip, fracs+1, derivs+1);
			if (*derivs) {
				return _deriv(arr, *fracs);
			} else {
				return _cubic_interp(arr, *fracs);
			};
		};
	};

};