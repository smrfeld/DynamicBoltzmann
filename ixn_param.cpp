#include "ixn_param.hpp"
#include "math.h"
#include "general.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Grid
	****************************************/

	/********************
	Constructor
	********************/

	Grid::Grid(std::string name, double min, double max, int n)
	{
		_name = name;
		_min = min;
		_max = max;
		_n = n;
		_delta = (max-min)/(n-1);
		_grid = new double[n];
		for (int i=0; i<n; i++) {
			_grid[i] = min+i*_delta;
		};
	};
	Grid::Grid(const Grid& other) 
	{
		_copy(other);
	};
	Grid& Grid::operator=(const Grid& other)
	{
		if (this != &other)
		{			
			_clean_up();
			_copy(other);
		};
		return *this;		
	};
	Grid::~Grid() {
		_clean_up();
	};
	void Grid::_copy(const Grid& other) {
		_name = other._name;
		_min = other._min;
		_max = other._max;
		_n = other._n;
		_delta = other._delta;
		_grid = new double[_n];
		std::copy( other._grid, other._grid + _n, _grid );
	};
	void Grid::_clean_up() {
		safeDelArr(_grid);
	};
	
	/********************
	Getters/setters
	********************/

	std::string Grid::name() const { return _name; };

	double Grid::delta() const { return _delta; };

	int Grid::n() const { return _n; };

	double Grid::get_by_idx(int i) const {
		return _grid[i];
	};

	/********************
	Surrounding idxs
	********************/

	int Grid::surrounding_idxs(double x) const
	{
		int i = (x - _min) / _delta;
		if (i==_n-1) {i--;};
		return i;
	};

	double Grid::frac_between(double x) const
	{
		return frac_between(x,surrounding_idxs(x));
	};

	double Grid::frac_between(double x, int i) const
	{
		return (x-_grid[i]) / _delta;
	};

	/********************
	Check if a given point is in the grid
	********************/

	bool Grid::in_grid(double x) const
	{
		if (x < _min || x > _max) { 
			return false; 
		} else {
			return true;
		};
	};
 
	/********************
	Print grid range
	********************/

	void Grid::print_grid_range() const
	{
		std::cout << "Grid: " << _name << " min: " << _min << " max: " << _max << std::endl;
	};

	/********************
	Write grid into an ofstream
	********************/

	void Grid::write_grid(std::string fname) const {
		std::ofstream f;
		f.open(fname);
		for (int i=0; i<_n; i++) {
			f << _grid[i] << "\n";
		};
		f.close();
	};

	/********************
	Test: if a sin func were defined on the grid
	********************/

	std::vector<double> Grid::test_sin() const {
		std::vector<double> x;
		for (int i=0; i<_n; i++) {
			x.push_back(sin(2*M_PI*_grid[i]/(_max-_min)));
		};
		return x;
	};
	std::vector<double> Grid::test_cos() const {
		std::vector<double> x;
		for (int i=0; i<_n; i++) {
			x.push_back(cos(2*M_PI*_grid[i]/(_max-_min)));
		};
		return x;	
	};

	/****************************************
	Ixn Param
	****************************************/

	/********************
	Constructor
	********************/

	IxnParam::IxnParam(std::string name, IxnParamType type, Species *sp, double min, double max, int n, double val0, int n_t) : IxnParam(name,type,sp,nullptr,min,max,n,val0,n_t) {};
	IxnParam::IxnParam(std::string name, IxnParamType type, Species *sp1, Species *sp2, double min, double max, int n, double val0, int n_t) : Grid(name,min,max,n)
	{
		_type = type;
		_sp1 = sp1;
		_sp2 = sp2;
		_val0 = val0;
		_n_t = n_t;

		_vals = new double[_n_t];
		std::fill_n(_vals,_n_t,0.0);
		_vals[0] = _val0;

		_asleep = new double[_n_t];
		std::fill_n(_asleep,_n_t,0.0);
		_awake = new double[_n_t];
		std::fill_n(_awake,_n_t,0.0);

		_bf = nullptr;
	};
	IxnParam::IxnParam(const IxnParam& other) : Grid(other) {
		_copy(other);
	};
	IxnParam& IxnParam::operator=(const IxnParam& other) {
		if (this != &other)
		{
			Grid::operator=(other);

			_clean_up();
			_copy(other);
		};
		return *this;
	};
	IxnParam::~IxnParam() {
		_clean_up();
	};
	void IxnParam::_copy(const IxnParam& other)
	{
		_type = other._type;
		_sp1 = other._sp1;
		_sp2 = other._sp2;
		_n_t = other._n_t;
		_vals = new double[_n_t];
		std::copy( other._vals, other._vals + _n_t, _vals );
		_val0 = other._val0;
		_asleep = new double[_n_t];
		std::copy( other._asleep, other._asleep + _n_t, _asleep );
		_awake = new double[_n_t];
		std::copy( other._awake, other._awake + _n_t, _awake );
		_bf = other._bf;
	};
	void IxnParam::_clean_up() {
		safeDelArr(_vals);
		safeDelArr(_asleep);
		safeDelArr(_awake);
	};


	/********************
	Set/check basis func pointer
	********************/

	void IxnParam::set_basis_func_ptr(BasisFunc* bf) {
		_bf = bf;
	};
	bool IxnParam::is_bf(BasisFunc *bf) {
		if (_bf==bf) { 
			return true;
		} else { 
			return false; 
		};
	};

	/********************
	Set IC
	********************/

	void IxnParam::set_init_cond(double val) {
		_val0 = val;
		_vals[0] = _val0;
	};

	/********************
	Validation
	********************/

	void IxnParam::validate_setup() const {
		std::cout << "--- Validate ixn param: " << name() << " ---" << std::endl; 
		if (_bf) {
			std::cout << "   Has basis func: " << _bf->name() << std::endl;
		} else {
			std::cerr << "ERROR: no basis func" << std::endl;
			exit(EXIT_FAILURE);
		};
	};

	/********************
	Getters/setters
	********************/

	double IxnParam::get_at_time(int it) const
	{
		return _vals[it];
	};

	/********************
	Calculate the next step
	********************/

	bool IxnParam::calculate_at_time(int it_next, double dt)
	{
		_vals[it_next] = _vals[it_next-1] + dt*_bf->get_at_time(it_next-1);
		return in_grid(_vals[it_next]);
	};

	/********************
	Moments from lattice
	********************/

	void IxnParam::moments_reset() 
	{
		for (int it=0; it<_n_t; it++) {
			_asleep[it] = 0.;
			_awake[it] = 0.;
		};
	};
	void IxnParam::moments_retrieve_at_time(MomentType moment_type, int it) {
		moments_retrieve_at_time(moment_type, it, 1);
	};
	void IxnParam::moments_retrieve_at_time(MomentType moment_type, int it, int batch_size)
	{
		if (_type == Hp) {
			if (_sp1) {
				if (moment_type==AWAKE) {
					_awake[it] += 1. * _sp1->count() / batch_size;
				} else if (moment_type==ASLEEP) {
					_asleep[it] += 1. * _sp1->count() / batch_size;
				};
			};
		} else if (_type == Jp) {
			if (_sp1 && _sp2) {
				if (moment_type==AWAKE) {
					_awake[it] += 1. * _sp1->nn_count(_sp2) / batch_size;
				} else if (moment_type==ASLEEP) {
					_asleep[it] += 1. * _sp1->nn_count(_sp2) / batch_size;
				};
			};
		};	
	};

	double IxnParam::moments_diff_at_time(int it) {
		return (_awake[it] - _asleep[it]);
	};

	/********************
	Write into an ofstream
	********************/

	void IxnParam::write_vals(std::string dir, int idx, int n_t_traj) const {
		std::ofstream f;
		f.open(dir+name()+"_"+pad_str(idx,4)+".txt");
		for (int i=0; i<n_t_traj; i++) {
			f << _vals[i];
			if (i != n_t_traj-1) { f << "\n"; };
		};
		f.close();
	};
	void IxnParam::write_vals(std::string dir, int idx1, int idx2, int n_t_traj) const {
		std::ofstream f;
		f.open(dir+name()+"_"+pad_str(idx1,4)+"_"+pad_str(idx2,2)+".txt");
		for (int i=0; i<n_t_traj; i++) {
			f << _vals[i];
			if (i != n_t_traj-1) { f << "\n"; };
		};
		f.close();
	};

	void IxnParam::write_moments(std::string dir, int idx, int n_t_traj) const {
		std::ofstream f;
		f.open(dir+name()+"_"+pad_str(idx,4)+".txt");
		for (int i=0; i<n_t_traj; i++) {
			f << _awake[i] << " " << _asleep[i];
			if (i != n_t_traj-1) { f << "\n"; };
		};
		f.close();
	};
	void IxnParam::write_moments(std::string dir, int idx1, int idx2, int n_t_traj) const {
		std::ofstream f;
		f.open(dir+name()+"_"+pad_str(idx1,4)+"_"+pad_str(idx2,2)+".txt");
		for (int i=0; i<n_t_traj; i++) {
			f << _awake[i] << " " << _asleep[i];
			if (i != n_t_traj-1) { f << "\n"; };
		};
		f.close();
	};

	/****************************************
	Array
	****************************************/

	Array::Array(std::vector<IxnParam*> ixn_params)
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
	Array::Array(IxnParam* ixn_param) : Array(std::vector<IxnParam*>{ixn_param}) {};
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

	void Array::write_vals(std::string dir, std::string name, int idx) const
	{
		std::ofstream f;
		f.open (dir+name+"_"+pad_str(idx,4)+".txt");
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
	Variational Term Trajectory
	****************************************/

	/********************
	Constructor
	********************/

	VarTerm::VarTerm(std::string name, IxnParam *num, BasisFunc *denom, std::vector<IxnParam*> denom_ixn_params, BasisFunc *num_bf, int n_ixn_params_in_num_bf, int n_t)
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
	VarTerm::VarTerm(const VarTerm& other)
	{
		_copy(other);
	};
	VarTerm& VarTerm::operator=(const VarTerm& other)
	{
		if (this != &other)
		{
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	VarTerm::~VarTerm()
	{
		_clean_up();
	};
	void VarTerm::_copy(const VarTerm& other) {
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
	void VarTerm::_clean_up() {
		safeDelArr(_num_bf_derivs);
	};

	/********************
	Set the pointers needed to update this term
	********************/

	void VarTerm::add_update_ptr(VarTerm* var_term)
	{
		_update_var_terms.push_back(var_term);
	};

	/********************
	Validation
	********************/

	void VarTerm::validate_setup() const {
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

	void VarTerm::calculate_at_time(int it_next, double dt)
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

	double VarTerm::get_at_time_by_idx(int it, int i) {
		return _vals[it].get_by_idx(i);
	};

	std::string VarTerm::name() {
		return _name;
	};

	/********************
	Write
	********************/

	void VarTerm::write_vals(std::string dir,int idx) const {
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

	/****************************************
	BasisFunc
	****************************************/

	/********************
	Constructor
	********************/

	BasisFunc::BasisFunc(std::string name, std::vector<IxnParam*> ixn_params) : Array(ixn_params) {
		_name = name;
		
		_derivs = new bool[_n_params];
		std::fill_n(_derivs,_n_params,false);

		_idxs_bounding = new int[_n_params];
		_idxs_p_cube = new int[_n_params];
		_idxs_ext_1 = new int[_n_params];
		_idxs_ext_2 = new int[_n_params];
		_fracs = new double[_n_params];
		_p_cube = new double[pow(4,_n_params)];

		_update_gathered = nullptr; // allocated later if needed
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
		_update_ptrs = other._update_ptrs;
		_derivs = new bool[_n_params];
		_idxs_bounding = new int[_n_params];
		_idxs_p_cube = new int[_n_params];
		_idxs_ext_1 = new int[_n_params];
		_idxs_ext_2 = new int[_n_params];
		_fracs = new double[_n_params];
		_p_cube = new double[pow(4,_n_params)];

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
	};

	/********************
	Add pointers needed to update
	********************/

	void BasisFunc::add_update_ptrs(IxnParam* ixn_param, VarTerm* var_term)
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

	void BasisFunc::update(int n_t, double dt, double dopt) 
	{
		// Go through all idxs
		for (int i=0; i<_val_len; i++) {
			// Go through all updating terms
			for (auto p: _update_ptrs) {
				// Go through all times
				for (int t=0; t<n_t; t++) {
					_vals[i] += dopt * dt * p.first->moments_diff_at_time(t) * p.second->get_at_time_by_idx(t, i);
				};
			};
		};
	};

	void BasisFunc::update_gather(int n_t, double dt, double dopt) 
	{
		if (!_update_gathered) {
			// alloc
			_update_gathered = new double[_val_len];
			std::fill_n(_update_gathered,_val_len,0.);
			std::cout << "!!! Alloced" << std::endl;
		};

		// Go through all idxs
		for (int i=0; i<_val_len; i++) {
			// Go through all updating terms
			for (auto p: _update_ptrs) {
				// Go through all times
				for (int t=0; t<n_t; t++) {
					_update_gathered[i] += dopt * dt * p.first->moments_diff_at_time(t) * p.second->get_at_time_by_idx(t, i);
				};
			};
		};
	};

	void BasisFunc::update_committ_gathered() 
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
			r *= exp(-pow(_ixn_params[ip]->get_at_time(it)-_ixn_params[ip]->get_by_idx(idxs[ip]),2)/(2.*1.0*_ixn_params[ip]->delta()))/sqrt(2.*M_PI*1.0*_ixn_params[ip]->delta());
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





