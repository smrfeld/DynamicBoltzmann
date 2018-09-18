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

	Domain1D::Domain1D(Iptr ixn_param, double min, double max, int no_pts) {
		_ixn_param = ixn_param;
		_min = min;
		_max = max;
		_no_pts = no_pts;
		_delta = (_max - _min) / (_no_pts - 1);
	};
	Domain1D::Domain1D(const Domain1D& other) {
		_copy(other);
	};
	Domain1D::Domain1D(Domain1D&& other) {
		_copy(other);
		other._reset();
	};
	Domain1D& Domain1D::operator=(const Domain1D& other) {
		if (this != &other)
		{
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	Domain1D& Domain1D::operator=(Domain1D&& other) {
		if (this != &other)
		{
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
		_min = other._min;
		_max = other._max;
		_delta = other._delta;
		_no_pts = other._no_pts;
	};
	void Domain1D::_reset()
	{
		_min = 0.0;
		_max = 0.0;
		_delta = 0.0;
		_no_pts = 0;
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
	int Domain1D::get_no_pts() const {
		return _no_pts;
	};
	double Domain1D::get_min() const {
		return _min;
	};
	double Domain1D::get_max() const {
		return _max;
	};
	double Domain1D::get_delta() const {
		return _delta;
	};

	// Get pt in domain
	double Domain1D::get_pt_by_idx(int i) const {
		if (i >= _no_pts) {
			std::cerr << ">>> Error: Domain1D::get_pt_by_idx <<< Idx: " << i << " is out of domain " << get_name() << " of size " << _no_pts << std::endl;
		};
		return _min + i * _delta;
	};

	// Check if point is in domain
	bool Domain1D::check_if_pt_is_inside_domain(double x) const {
		if (x < _min || x > _max) { 
			return false; 
		} else {
			return true;
		};
	};

	// Get indexes surrounding a point
	// ie point is between i and i+1 where i is returned
	int Domain1D::get_idxs_surrounding_pt(double x) const {
		int i = (x - _min) / _delta;
		if (i==_no_pts-1) {i--;};
		return i;
	};

	// Get fraction of a point between successive points
	double Domain1D::get_frac_between(double x) const {
		return get_frac_between(x,get_idxs_surrounding_pt(x));
	};
	// Second optional specification: the return of the surrounding idxs
	double Domain1D::get_frac_between(double x, int i) const {
		return (x - get_pt_by_idx(i)) / _delta;
	};

	// Print domain range
	void Domain1D::print_bounds() const {
		std::cout << "Domain: " << get_name() << " min: " << _min << " max: " << _max << std::endl;
	};





























	/****************************************
	Array
	****************************************/

	Array::Array(std::vector<Domain1D> domain)
	{
		_domain = domain;
		_n_params = _domain.size();

		// Values
		_val_len = 1;
		for (auto v: _domain) { _val_len *= v.get_no_pts(); };
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
				pwr *= _domain[jv].get_no_pts();
			};
			_dim_pwrs.push_back(pwr); 
		};
	};
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
		_domain = other._domain;
		_dim_pwrs = other._dim_pwrs;
		_vals = new double[_val_len];
		std::copy( other._vals, other._vals + _val_len, _vals );
		// reset other
		other._n_params = 0;
		other._val_len = 0;
		other._domain.clear();
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
			_domain = other._domain;
			_dim_pwrs = other._dim_pwrs;
			_vals = new double[_val_len];
			std::copy( other._vals, other._vals + _val_len, _vals );
			// reset other
			other._n_params = 0;
			other._val_len = 0;
			other._domain.clear();
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
		_domain = other._domain;
		_dim_pwrs = other._dim_pwrs;

		_vals = new double[_val_len];
		std::copy( other._vals, other._vals + _val_len, _vals );
	};
	void Array::_clean_up() {
		safeDelArr(_vals);
	};

	/********************
	Get domain
	********************/

	const std::vector<Domain1D>& Array::get_domain() const {
		return _domain;
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
				f << _domain[ip].get_pt_by_idx(idxs[ip]);
				if (ip != _n_params-1) { f << " "; };
			};
			f << "\n";
		};

		f.close();

		// Clean up
		safeDelArr(idxs);
	};

	void Array::write_vals(std::string fname) const
	{
		std::ofstream f;
		f.open (fname);
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
			if (_domain[id].get_name() != other._domain[id].get_name()) {
				return false;
			};
		};
		return true;
	};

	/********************
	Zero
	********************/

	void Array::reset_to_zero()
	{
		std::fill_n(_vals,_val_len,0.);
	};









































	/****************************************
	DiffEqRHS
	****************************************/

	/********************
	Constructor
	********************/

	DiffEqRHS::DiffEqRHS(std::string name, std::vector<Domain1D> domain) : Array(domain) {
		_name = name;

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
	DiffEqRHS::DiffEqRHS(const DiffEqRHS& other) : Array(other) {
		_copy(other);
	};
	DiffEqRHS& DiffEqRHS::operator=(const DiffEqRHS& other) {
		if (this != &other)
		{
			Array::operator=(other);

			_clean_up();

			_copy(other);
		};
		return *this;
	};
	DiffEqRHS::~DiffEqRHS() {
		_clean_up();
	};
	void DiffEqRHS::_copy(const DiffEqRHS& other) {
		_name = other._name;

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
	void DiffEqRHS::_clean_up() {
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

	/********************
	Set prev nesterov
	********************/

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

	/********************
	Validate
	********************/

	// Validate setup
	void DiffEqRHS::validate_setup() const {
		std::cout << "--- Validate Basis func " << _name << " ---" << std::endl;
	};

	/********************
	Get values, if they are in the lattice
	********************/

	// Recursion to fill p
	// idx_deltas will in the inner most loop be -1, 0, 1, 2 for i-1,i,i+1,i+2
	// Initial dim parameter = 0
	// idxs = contains the i
	void DiffEqRHS::_fill_p(int dim)
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
					} else if (idx > _domain[d].get_no_pts()-1) {
						outside = true; // Yes it's outside on at least one dim
						_idxs_ext_1[d] = _domain[d].get_no_pts()-1;
						_idxs_ext_2[d] = _domain[d].get_no_pts()-2;
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
	void DiffEqRHS::_get_bounding(int it, bool safe)
	{
		// Check that x is in the box
		if (safe) {
			for (int id=0; id<_n_params; id++) {
				if (!(_domain[id].check_if_pt_is_inside_domain(_domain[id].get_ixn_param()->get_val_at_timepoint(it)))) {
					std::cerr << "ERROR - " << _domain[id].get_ixn_param()->get_val_at_timepoint(it) << " is outside the grid:" << std::endl;
					_domain[id].print_bounds();
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
			x = _domain[id].get_ixn_param()->get_val_at_timepoint(it);
			_idxs_bounding[id] = _domain[id].get_idxs_surrounding_pt(x);
			_fracs[id] = _domain[id].get_frac_between(x,_idxs_bounding[id]);
		};

		// Get bounding box
		_fill_p(0);
	};

	/********************
	Get values, if they are in the lattice
	********************/

	double DiffEqRHS::get_val_at_timepoint(int it) 
	{
		// Create the bounding box
		_get_bounding(it);

		return _n_cubic_interp(_n_params,_p_cube,_fracs,_derivs);
	};

	double DiffEqRHS::get_deriv_at_timepoint(int it, int i_dim) {
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

	std::string DiffEqRHS::get_name() const {
		return _name;
	};

	/********************
	Calculate the new basis function
	********************/

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
			};
		};
	};

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

	/****************************************
	DiffEqRHS - PRIVATE
	****************************************/

	/********************
	Interpolation
	********************/

	double DiffEqRHS::_cubic_interp(double p[4], double f) {
		return p[1] + 0.5 * f *(p[2] - p[0] + f*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + f*(3.0*(p[1] - p[2]) + p[3] - p[0])));
	};
	double DiffEqRHS::_deriv(double p[4], double f) {
		return 0.5 * (p[2] - p[0]) + 2.0 * f * (p[0] - 2.5 * p[1] + 2.0 * p[2] - 0.5 * p[3]) + 3.0 * f * f *( - 0.5 * p[0] + 1.5 * p[1] - 1.5 * p[2] + 0.5 * p[3]);
	};
	double DiffEqRHS::_n_cubic_interp(int dim, double* p, double fracs[], bool derivs[])
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