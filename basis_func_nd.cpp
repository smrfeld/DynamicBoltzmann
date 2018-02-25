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

	Dim::Dim(std::string nameIn, double minIn, double maxIn, int nIn) {
		n = nIn;
		min = minIn;
		max = maxIn;
		name = nameIn;
		delta = (max-min)/(n-1);
		grid = new double[n];
		for (int i=0; i<n; i++) {
			this->grid[i] = min + i*delta;
		};
	};
	Dim::Dim(const Dim& d) 
	{
		n = d.n;
		min = d.min;
		max = d.max;
		delta = d.delta;
		name = d.name;
		grid = new double[n];
		std::copy( d.grid, d.grid + n, grid );
	};
	Dim& Dim::operator=(const Dim& d)
	{
		if (this != &d)
		{
			n = d.n;
			min = d.min;
			max = d.max;
			delta = d.delta;
			name = d.name;

			safeDelArr(grid);

			grid = new double[n];
			std::copy( d.grid, d.grid + n, grid );
		};
		return *this;		
	};
	Dim::~Dim() {
		safeDelArr(grid);
	};

	/****************************************
	GridND
	****************************************/
	GridND::GridND(int n_dim, std::vector<Dim> dims)
	{
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
	void GridND::_copy(const GridND& g) {
		_n_dim = g._n_dim;
		_val_len = g._val_len;
		_dims = g._dims;

		_vals = new double[_val_len];
		std::copy( g._vals, g._vals + _val_len, _vals );
	};

	/********************
	Get/set an element by index
	********************/

	double GridND::get(int *idxs)
	{
		int i=0;
		for (int id=0; id<_n_dim; id++) {
			i += _dim_pwrs[id] * idxs[id];
		};
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

	void GridND::write_to_file(std::string fname) 
	{
		std::ofstream f;
		f.open (fname);

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

	/****************************************
	BasisFuncND
	****************************************/

	/********************
	Constructor
	********************/

	BasisFuncND::BasisFuncND(int n_dim, std::vector<Dim> dims) : GridND(n_dim,dims) {
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

	double BasisFuncND::get(int *idxs) {
		return GridND::get(idxs);
	};

	void BasisFuncND::set(int *idxs, double val) {
		GridND::set(idxs,val);
	};

	void BasisFuncND::write_to_file(std::string fname) {
		GridND::write_to_file(fname);
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


