#include "basis_func_nd.hpp"
#include "math.h"
#include "general.hpp"
#include <algorithm>
#include <iostream>
#include <fstream>

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Dimension information
	****************************************/

	Dim::Dim(std::string nameIn, double minIn, double maxIn, int nIn) : Dim(nameIn, minIn, maxIn, nIn, NONE, "", "") {};
	Dim::Dim(std::string nameIn, double minIn, double maxIn, int nIn, DimType typeIn, std::string speciesIn) : Dim(nameIn, minIn, maxIn, nIn, typeIn, speciesIn, "") {};
	Dim::Dim(std::string nameIn, double minIn, double maxIn, int nIn, DimType typeIn, std::string species1In, std::string species2In) 
	{
		name = nameIn;
		
		type = typeIn;

		n = nIn;
		min = minIn;
		max = maxIn;
		delta = (max-min)/(n-1);
		grid = new double[n];
		for (int i=0; i<n; i++) {
			this->grid[i] = min + i*delta;
		};
		
		species1 = species1In;
		species2 = species2In;
		sp1 = nullptr;
		sp2 = nullptr;

		awake = 0.0;
		asleep = 0.0;

		// Default: use all dimensions as basis
		bf_use_all_dim = true;
	};
	Dim::Dim(const Dim& d) 
	{
		_copy(d);
	};
	Dim& Dim::operator=(const Dim& d)
	{
		if (this != &d)
		{			
			safeDelArr(grid);
			_copy(d);
		};
		return *this;		
	};
	Dim::~Dim() {
		safeDelArr(grid);
	};

	/********************
	Copy function
	********************/

	void Dim::_copy(const Dim& d) {
		name = d.name;
		
		type = d.type;

		n = d.n;
		min = d.min;
		max = d.max;
		delta = d.delta;

		grid = new double[n];
		std::copy( d.grid, d.grid + n, grid );

		species1 = d.species1;
		species2 = d.species2;
		if (d.sp1) {
			sp1 = d.sp1;
		};
		if (d.sp2) {
			sp2 = d.sp2;
		};

		awake = d.awake;
		asleep = d.asleep;

		bf_use_all_dim = d.bf_use_all_dim;
		bf_n_dim = d.bf_n_dim;
		bf_dim_idxs = d.bf_dim_idxs;
		bf_dim_names = d.bf_dim_names;
	};

	/********************
	Specify dimensions to use as basis
	********************/

	void Dim::set_basis_func_dims(std::string dim_name) 
	{
		bf_use_all_dim = false;
		auto it = std::find(bf_dim_names.begin(), bf_dim_names.end(), dim_name);
		if (it == bf_dim_names.end()) {
			bf_dim_names.push_back(dim_name);
			bf_n_dim = bf_dim_names.size();
		};
	};
	void Dim::set_basis_func_dims(std::vector<std::string> dim_name) {
		for (auto s: dim_name) {
			set_basis_func_dims(s);
		};
	};

	/********************
	Update moments from species
	********************/

	void Dim::append_moments_from_latt(MomentType moment_type) {
		if (type == H) {
			if (sp1) {
				if (moment_type==AWAKE) {
					awake += sp1->count;
				} else if (moment_type==ASLEEP) {
					asleep += sp1->count;
				};
			};
		} else if (type == J) {
			if (sp1 && sp2) {
				if (moment_type==AWAKE) {
					awake += sp1->nn_count[sp2];
				} else if (moment_type==ASLEEP) {
					asleep += sp1->nn_count[sp2];
				};
			};
		};
	};

	/********************
	Comparators
	********************/

	bool Dim::operator==(const Dim& d) 
	{
		return std::tie(name,n,min,max) == std::tie(d.name,d.n,d.min,d.max);
	};
	bool Dim::operator!=(const Dim& d) {
		return !(*this==d);
	}

	/****************************************
	GridND
	****************************************/

	GridND::GridND(std::string name, int n_dim, std::vector<Dim> dims)
	{
		this->name = name;

		_n_dim = n_dim;
		_dims = dims;

		// Values
		_val_len = 1;
		for (auto v: _dims) { _val_len *= v.n; };
		_vals = new double[_val_len];

		// Zero by default
		std::fill_n(_vals,_val_len,0.0);

		// Update dimension powers
		int pwr;
		for (int iv=0; iv<_n_dim; iv++) 
		{
			pwr = 1;
			for (int jv=iv+1; jv<_n_dim; jv++)
			{
				pwr *= _dims[jv].n;
			};
			_dim_pwrs.push_back(pwr); 
		};
	};
	GridND::GridND(std::string name, Dim dim) : GridND(name,1,std::vector<Dim>{dim}) {};
	GridND::GridND(const GridND& g) 
	{
		_copy(g);
	};
	GridND& GridND::operator=(const GridND& g)
	{
		if (this != &g) {
			safeDelArr(_vals);
			_copy(g);	
		};
		return *this;
	};
	GridND::~GridND()
	{
		safeDelArr(_vals);
	};

	/********************
	Internal copy
	********************/

	void GridND::_copy(const GridND& g) {
		_n_dim = g._n_dim;
		_val_len = g._val_len;
		_dims = g._dims;
		_dim_pwrs = g._dim_pwrs;
		this->name = g.name;

		_vals = new double[_val_len];
		std::copy( g._vals, g._vals + _val_len, _vals );
	};

	/********************
	Get/set an element by index
	********************/

	double GridND::get(int *idxs) const
	{
		int i=0;
		for (int id=0; id<_n_dim; id++) {
			i += _dim_pwrs[id] * idxs[id];
		};
		return _vals[i];
	};

	double GridND::get(int i) const
	{
		return _vals[i];
	};

	void GridND::set(int *idxs, double val)
	{
		int i=0;
		for (int id=0; id<_n_dim; id++) {
			i += _dim_pwrs[id] * idxs[id];
		};
		_vals[i] = val;
	};

	void GridND::set(int i, double val)
	{
		_vals[i] = val;
	};

	/********************
	Multiply all values by some constant
	********************/

	void GridND::multiply_all(double val)
	{
		for (int i=0; i<_val_len; i++) {
			_vals[i] *= val;
		};
	};

	/********************
	Increment by a whole grid
	********************/

	void GridND::increment(const GridND& g)
	{
		/*
		if (!check_dims(g)) { 
			std::cerr << "Could not increment because dims don't match" << std::endl;
			exit(EXIT_FAILURE); 
		};
		*/

		for (int i=0; i<_val_len; i++) {
			_vals[i] += g.get(i);
		};
	};

	/********************
	Get indexes by element
	********************/

	int* GridND::get_idxs(int i) {
		int *idxs = new int[_n_dim];
		int rem = i;
		for (int id=0; id<_n_dim; id++) {
			idxs[id] = int(rem/_dim_pwrs[id]);
			rem = rem % _dim_pwrs[id];
		};
		return idxs;
	};

	/********************
	Write to file
	********************/

	void GridND::write_to_file(std::string dir, int idx) 
	{
		std::ofstream f;
		f.open (dir + name + "_" + pad_str(idx,4));

		// Go through everything to write
		int *idxs = new int[_n_dim];
		for (int i=0; i<_val_len; i++) {
			idxs = get_idxs(i);
			for (int j=0; j<_n_dim; j++) {
				f << _dims[j].grid[idxs[j]] << " ";
			};
			f << " " << _vals[i] << "\n";
		};
		f.close();
		safeDelArr(idxs);
	};

	/********************
	Check dimensions against another grid
	********************/

	bool GridND::check_dims(const GridND& g)
	{
		// Check no dims
		if (_n_dim != g._n_dim) { 
			std::cerr << "ERROR! No dimensions don't match in update step" << std::endl;
			return false;
		};

		// Check each dim
		for (int id=0; id<_n_dim; id++) {
			if (_dims[id] != g._dims[id]) {
				return false;
			};
		};
		return true;
	};

	/********************
	Zero
	********************/

	void GridND::zero()
	{
		std::fill_n(_vals,_val_len,0.);
	};

	/****************************************
	BasisFuncND
	****************************************/

	/********************
	Constructor
	********************/

	BasisFuncND::BasisFuncND(std::string name, int n_dim, std::vector<Dim> dims) : GridND(name,n_dim,dims) {
		_idxs_bounding = new int[_n_dim];
		_idxs_p_cube = new int[_n_dim];
		_idxs_ext_1 = new int[_n_dim];
		_idxs_ext_2 = new int[_n_dim];
		_fracs = new double[_n_dim];
		_p_cube = new double[pow(4,_n_dim)];
	};
	BasisFuncND::BasisFuncND(const BasisFuncND& bf) : GridND(bf) {
		_copy(bf);
	};
	BasisFuncND& BasisFuncND::operator=(const BasisFuncND& bf) {
		if (this != &bf)
		{
			GridND::operator=(bf);

			_clean_up();

			_copy(bf);
		};
		return *this;
	};
	BasisFuncND::~BasisFuncND() {
		_clean_up();
	};
	void BasisFuncND::_copy(const BasisFuncND& bf) {
		_idxs_bounding = new int[_n_dim];
		_idxs_p_cube = new int[_n_dim];
		_idxs_ext_1 = new int[_n_dim];
		_idxs_ext_2 = new int[_n_dim];
		_fracs = new double[_n_dim];
		_p_cube = new double[pow(4,_n_dim)];

		std::copy( bf._idxs_bounding, bf._idxs_bounding + _n_dim, _idxs_bounding );
		std::copy( bf._idxs_p_cube, bf._idxs_p_cube + _n_dim, _idxs_p_cube );
		std::copy( bf._idxs_ext_1, bf._idxs_ext_1 + _n_dim, _idxs_ext_1 );
		std::copy( bf._idxs_ext_2, bf._idxs_ext_2 + _n_dim, _idxs_ext_2 );
		std::copy( bf._fracs, bf._fracs + _n_dim, _fracs );
		std::copy( bf._p_cube, bf._p_cube + int(pow(4,_n_dim)), _p_cube );
	};
	void BasisFuncND::_clean_up() {
		safeDelArr(_idxs_bounding);
		safeDelArr(_idxs_p_cube);
		safeDelArr(_idxs_ext_1);
		safeDelArr(_idxs_ext_2);
		safeDelArr(_fracs);
		safeDelArr(_p_cube);
	};

	/********************
	Get values, if they are in the lattice
	********************/

	// Recursion to fill p
	// idx_deltas will in the inner most loop be -1, 0, 1, 2 for i-1,i,i+1,i+2
	// Initial dim parameter = 0
	// idxs = contains the i
	void BasisFuncND::_fill_p(int dim)
	{
		// i-1 to i+2
		for(_idxs_p_cube[dim] = 0; _idxs_p_cube[dim] < 4; ++_idxs_p_cube[dim]) 
		{
			if (dim != _n_dim-1)
			{
				// Inception - need another for loop
				_fill_p(dim+1);
			} else {
				// Do something!

				// Print
				/*
				std::cout << "p cube idxs: ";
				for (int a=0; a<_n_dim; a++) {
					std::cout << _idxs_p_cube[a] << " ";
				};
				*/

				// Check if point is outside box
				bool outside=false; // Is it outside on at least one dim?
				int idx;
				for (int d=0; d<_n_dim; d++) {
					idx = _idxs_bounding[d]+_idxs_p_cube[d]-1;
					if (idx < 0) {
						outside = true; // Yes it's outside on at least one dim
						_idxs_ext_1[d] = 0;
						_idxs_ext_2[d] = 1;
					} else if (idx > _dims[d].n-1) {
						outside = true; // Yes it's outside on at least one dim
						_idxs_ext_1[d] = _dims[d].n-1;
						_idxs_ext_2[d] = _dims[d].n-2;
					} else {
						_idxs_ext_1[d] = idx;
						_idxs_ext_2[d] = idx;
					};
				};

				// Index in _p_cube to write to
				int i_write=0;
				for (int d=0; d<_n_dim; d++) {
					i_write += pow(4,_n_dim-d-1) * _idxs_p_cube[d];
				};

				// Finally, write
				if (outside) {
					_p_cube[i_write] = 2*get(_idxs_ext_1) - get(_idxs_ext_2);
					// std::cout << "writing: 2*(" << _idxs_ext_1[0] << "," << _idxs_ext_1[1] << ") - (" << _idxs_ext_2[0] << "," << _idxs_ext_2[1] << ")";
				} else {
					_p_cube[i_write] = get(_idxs_ext_1);
					// std::cout << "writing: (" << _idxs_ext_1[0] << "," << _idxs_ext_1[1] << ")";
				};
				// std::cout << std::endl;
			};
		};
	};

	// Function to get bounding cube of 4 points
	void BasisFuncND::_get_bounding(double *x) 
	{
		// Check that x is in the box
		for (int id=0; id<_n_dim; id++) {
			if (x[id] < _dims[id].min || x[id] > _dims[id].max) {
				std::cerr << "ERROR - " << x[id] << " is outside the interval " << _dims[id].min << " to " << _dims[id].max << std::endl;
				exit(EXIT_FAILURE);
			};
		};

		// Get the interval - i-1,i,POINT,i+1,i+2
		// indexes in a dimension of length n run from 0 to n-1
		// this should return i = 0 to n-2
		// left boundary: i=0 -> -1,0,1,2 where point i=-1 needs extrapolation later
		// right boundary: i=n-2 -> n-3,n-2,n-1,n where point i=n needs extrapolation later
		// ALSO: get the fraction this point is between two successive, i.e. between i and i+1
		for (int id=0; id<_n_dim; id++) {
			_idxs_bounding[id] = (x[id] - _dims[id].min) / _dims[id].delta;
			if (_idxs_bounding[id] == _dims[id].n-1) { _idxs_bounding[id]--; };
			_fracs[id] = (x[id] - _dims[id].grid[_idxs_bounding[id]]) / _dims[id].delta;
		};

		// Get bounding box
		_fill_p(0);
	};

	/********************
	Get values, if they are in the lattice
	********************/

	double BasisFuncND::get_val(double *x) 
	{
		// No derivs
		bool derivs[_n_dim];
		std::fill_n(derivs,_n_dim,false);

		return get_deriv(x,derivs);
	};

	double BasisFuncND::get_deriv(double *x, bool *derivs) {
		// Create the bounding box
		_get_bounding(x);

		// Interpolate
		return _n_cubic_interp(_n_dim,_p_cube,_fracs,derivs);
	};

	/********************
	Update with delta f
	********************/

	void BasisFuncND::update(const GridND &g) 
	{
		// Check dimensions
		if (!check_dims(g)) { 
			std::cerr << "ERROR! Could not update BF b/c dims don't match" << std::endl; 
			exit(EXIT_FAILURE); 
		};

		// Update
		for (int i=0; i<_val_len; i++) {
			_vals[i] += g.get(i);
		};
	};

	/********************
	Test fill in various dimensions
	********************/

	double BasisFuncND::test_func_2d(double x, double y)
	{
		return sin(2*M_PI*x/(_dims[0].max-_dims[0].min))*cos(2*M_PI*y/(_dims[1].max-_dims[1].min));
	};
	void BasisFuncND::test_fill_2d() {
		if (_n_dim != 2) { return; };

		int idxs[2];
		for (int i=0; i<_dims[0].n; i++) {
			for (int j=0; j<_dims[1].n; j++) {
				idxs[0] = i;
				idxs[1] = j;
				set(idxs, test_func_2d(_dims[0].grid[i],_dims[1].grid[j]));
			};
		};
	};
	double BasisFuncND::test_func_3d(double x, double y, double z)
	{
		return sin(2*M_PI*x/(_dims[0].max-_dims[0].min))*cos(2*M_PI*y/(_dims[1].max-_dims[1].min)) + pow(z,2)*y;
	};
	void BasisFuncND::test_fill_3d() {
		if (_n_dim != 3) { return; };

		int idxs[3];
		for (int i=0; i<_dims[0].n; i++) {
			for (int j=0; j<_dims[1].n; j++) {
				for (int k=0; k<_dims[2].n; k++) {
					idxs[0] = i;
					idxs[1] = j;
					idxs[2] = k;
					set(idxs, test_func_3d(_dims[0].grid[i],_dims[1].grid[j],_dims[2].grid[k]));
				};
			};
		};
	};

	/********************
	From parent
	********************/

	double BasisFuncND::get(int *idxs) const {
		return GridND::get(idxs);
	};

	double BasisFuncND::get(int i) const {
		return GridND::get(i);
	};

	void BasisFuncND::set(int *idxs, double val) {
		GridND::set(idxs,val);
	};

	void BasisFuncND::write_to_file(std::string dir, int idx) {
		GridND::write_to_file(dir, idx);
	}; 

	/****************************************
	BasisFuncND - PRIVATE
	****************************************/

	/********************
	Interpolation
	********************/

	double BasisFuncND::_cubic_interp(double p[4], double f) {
		return p[1] + 0.5 * f *(p[2] - p[0] + f*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + f*(3.0*(p[1] - p[2]) + p[3] - p[0])));
	};
	double BasisFuncND::_deriv(double p[4], double f) {
		return 0.5 * (p[2] - p[0]) + 2.0 * f * (p[0] - 2.5 * p[1] + 2.0 * p[2] - 0.5 * p[3]) + 3.0 * f * f *( - 0.5 * p[0] + 1.5 * p[1] - 1.5 * p[2] + 0.5 * p[3]);
	};
	double BasisFuncND::_n_cubic_interp(int dim, double* p, double fracs[], bool derivs[])
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


