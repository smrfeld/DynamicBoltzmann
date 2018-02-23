#include "basis_func.hpp"
#include <iostream>
#include <algorithm>
#include "general.hpp"
#include <fstream>
#include "math.h"

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	BasisFunc1D
	****************************************/

	// Constructor
	BasisFunc1D::BasisFunc1D(double nu_min, double nu_max, int n_nu)
	{
		// Check min no pts is satisfied
		if (n_nu <= 4) { 
			std::cerr << "ERROR - n_nu > 4 required." << std::endl; 
			exit(EXIT_FAILURE); 
		};

		this->_min = nu_min;
		this->_max = nu_max;
		this->_n = n_nu;
		this->_delta = (nu_max - nu_min) / (n_nu - 1);

		// All zeros by default
		this->_vals = new double[n_nu];
		std::fill_n(this->_vals, n_nu, 0.0);
		/*
		for (int i=0; i<_n; i++) {
			_vals[i] = 0.1*sin(2.0*M_PI*i/(_n-1))+0.1;
		};
		*/
		this->_deriv = new double[n_nu];
		// Update derivs
		_update_derivs();

		// TEMP: fill the basis func with sin
		// test_fill_sin();
	};
	BasisFunc1D::BasisFunc1D(const BasisFunc1D& bf)
	{
		// Copy standard
		this->_min = bf._min;
		this->_max = bf._max;
		this->_n = bf._n;
		this->_delta = bf._delta;

	    // Allocate then copy soln, nus
	    this->_vals = new double[bf._n];
	    *(this->_vals) = *(bf._vals);
	    this->_deriv = new double[bf._n];
	    *(this->_deriv) = *(bf._deriv);
	};
	BasisFunc1D& BasisFunc1D::operator=(const BasisFunc1D& bf)
	{
		if (this != &bf)  
		{
			// Copy standard
			this->_min = bf._min;
			this->_max = bf._max;
			this->_n = bf._n;
			this->_delta = bf._delta;

			// Free the existing resource.  
			safeDelArr(this->_vals);  
			safeDelArr(this->_deriv);  

		    // Allocate then copy soln, nus
		    this->_vals = new double[bf._n];
		    *(this->_vals) = *(bf._vals);
		    this->_deriv = new double[bf._n];
		    *(this->_deriv) = *(bf._deriv);
		};
		return *this;
	};
	BasisFunc1D::~BasisFunc1D() {
		safeDelArr(this->_vals);
		safeDelArr(this->_deriv);  
	};

	/********************
	Get values
	********************/

	double BasisFunc1D::get_val(double nu) 
	{
		return _cubic_interp(_vals, nu);
	};
	double BasisFunc1D::get_deriv(double nu) {
		return _cubic_interp(_deriv, nu);
	};

	/********************
	Write lattice to a file
	********************/
	
	void BasisFunc1D::write_to_file(std::string fname) {
		_write_to_file(_vals,fname);
	};
	void BasisFunc1D::write_to_file(std::string fname, int npts) {
		_write_to_file(_vals,fname,npts);
	};
	void BasisFunc1D::write_deriv_to_file(std::string fname) {
		_write_to_file(_deriv,fname);
	};
	void BasisFunc1D::write_deriv_to_file(std::string fname, int npts) {
		_write_to_file(_deriv,fname,npts);
	};

	/********************
	IO
	********************/

	std::ostream& operator<<(std::ostream& os,  BasisFunc1D& bf)
	{
		for (int i=0; i<bf._n; i++) 
		{ 
			os << bf._min + i*bf._delta << " " << bf._vals[i];
			if (i != bf._n-1) {
				os << "\n";
			};
		};
		return os;
	};

	/********************
	Accessor
	********************/

	double& BasisFunc1D::operator[](std::size_t idx) { return _vals[idx]; };

	/********************
	Update the values
	********************/

	void BasisFunc1D::update(double *dvals) {
		for (int i=0; i<_n; i++) {
			_vals[i] += dvals[i];
		};

		// Update the deriv
		_update_derivs();
	};

	/********************
	Update the values, smoothing delta F before
	********************/
	void BasisFunc1D::update_smoothed(double *dvals)
	{
		// Smooth the updates by an average
		// This helps to get rid of edge effects
		int window = 11;
		double *sm = new double[_n];
		// Left edge
		for (int i=0; i<int(window/2); i++) {
			sm[i] = 0.;
			for (int j=0; j<int(window/2)+1+i; j++) {
				sm[i] += dvals[j];
			};
			sm[i] /= int(window/2)+1+i;
		};
		// Interior
		for (int i=int(window/2); i<_n-int(window/2); i++) {
			sm[i] = 0.;
			for (int j=i-int(window/2); j<=i+int(window/2); j++) {
				sm[i] += dvals[j];
			}
			sm[i] /= window;
		};
		// Right edge
		for (int i=0; i<int(window/2); i++) {
			sm[_n-1-i] = 0.;
			for (int j=0; j<int(window/2)+1+i; j++) {
				sm[_n-1-i] += dvals[_n-1-j];
			};
			sm[_n-1-i] /= int(window/2)+1+i;
		};

		// Update all the values
		for (int i=0; i<_n; i++) {
			_vals[i] += sm[i];
		};

		// Clear
		delete[] sm;


		// Update the deriv
		_update_derivs();
	};

	/********************
	Test: fill with a sin func
	********************/

	void BasisFunc1D::test_fill_sin() {
	    for (int i=0; i<_n; i++) {
			this->_vals[i] = sin(2.0*M_PI*i/(_n-1))+0.1;
		};
		_update_derivs();
	};

	/****************************************
	BasisFunc1D - PRIVATE
	****************************************/

	// For interpolation, see:
	// http://www.paulinternet.nl/?page=bicubic

	/********************
	Get the two points left and right of a given point
	********************/

	double *BasisFunc1D::_get_interval(double *s, double x) {
		if (x < _min || x > _max) { 
			std::cerr << "ERROR - outside interval." << std::endl; 
			exit(EXIT_FAILURE); 
		};		
		int i = (x - _min)/_delta;
		if (i == _n - 1) { i--; };
		static double p[4];
		p[1] = s[i];
		p[2] = s[i+1];

		if (i == 0) {
			// Between the first two points - make line for first point
			p[0] = 2.*s[i] - s[i+1];
			p[3] = s[i+2];
		} else if (i == _n - 2) {
			// Between last two points - make line for last point
			p[0] = s[i-1];
			p[3] = 2.*s[i+1]-s[i];
		} else {
			p[0] = s[i-1];
			p[3] = s[i+2];
		};
		return p;
	};

	// Note: the corresponding nu values to p0,p1,p2,p3
	// are x=-1,0,1,2 so the interp value is converted to 0-1
	double BasisFunc1D::_cubic_interp(double *s, double x)
	{
		double *p = _get_interval(s, x);
		double f = fmod(x,_delta)/_delta;
		if (f > 1.0-1E-6) { f=0.0; }; // always use 0 rather than 1
		if (_max - x < 1E-6) { f=1.0; }; // very end
		return p[1] + 0.5 * f*(p[2] - p[0] + f*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + f*(3.0*(p[1] - p[2]) + p[3] - p[0])));
	};

	double BasisFunc1D::_get_deriv(double x)
	{
		double *p = _get_interval(_vals,x);
		double f = fmod(x,_delta)/_delta;
		if (f > 1.0-1E-6) { f=0.0; }; // always use 0 rather than 1
		if (_max - x < 1E-6) { f=1.0; }; // very end
		return 0.5 * (p[2] - p[0]) + 2.0 * f * (p[0] - 2.5 * p[1] + 2.0 * p[2] - 0.5 * p[3]) + 3.0 * f * f *( - 0.5 * p[0] + 1.5 * p[1] - 1.5 * p[2] + 0.5 * p[3]);
	};

	/********************
	Write to a file
	********************/

	void BasisFunc1D::_write_to_file(double *s, std::string fname) 
	{
		std::ofstream f;
		f.open (fname);
		for (int i=0; i<_n; i++) {
			f << _min+i*_delta << " " << s[i];
			if (i != _n-1) { f << "\n";};
		};
		f.close();
	};

	void BasisFunc1D::_write_to_file(double *s, std::string fname, int npts)
	{
		std::ofstream f;
		f.open (fname);
		double d = (_max-_min) / (npts-1);
		for (int i=0; i<npts; i++) {
			f << _min + i*d << " " << _cubic_interp(s, _min + i*d);
			if (i != npts-1) { f << "\n";};
		};
		f.close();
	};

	/********************
	Calculate derivative at all points and store
	********************/

	void BasisFunc1D::_update_derivs() {
		for (int i=0; i<_n; i++) {
			_deriv[i] = _get_deriv(_min+i*_delta);
		};
	};















	/****************************************
	BasisFunc2D
	****************************************/

	// Constructor
	BasisFunc2D::BasisFunc2D(double nu_min1, double nu_max1, int n_nu1, double nu_min2, double nu_max2, int n_nu2)
	{
		// Check min no pts is satisfied
		if (n_nu1 <= 4 || n_nu2 <= 4) { 
			std::cerr << "ERROR - n_nu > 4 required." << std::endl; 
			exit(EXIT_FAILURE); 
		};

		this->_min1 = nu_min1;
		this->_max1 = nu_max1;
		this->_n1 = n_nu1;
		this->_delta1 = (nu_max1 - nu_min1) / (n_nu1 - 1);
		this->_min2 = nu_min2;
		this->_max2 = nu_max2;
		this->_n2 = n_nu2;
		this->_delta2 = (nu_max2 - nu_min2) / (n_nu2 - 1);

		// All zeros by default
		this->_vals = new double*[_n1];
		for (int i=0; i<_n1; i++) {
			_vals[i] = new double[_n2];
			std::fill_n(_vals[i], _n2, 0.0);
		};
		this->_valsb = new double*[_n1+2];
		for (int i=0; i<_n1+2; i++) {
			_valsb[i] = new double[_n2+2];
			std::fill_n(_valsb[i], _n2+2, 0.0);
		};
	};
	BasisFunc2D::BasisFunc2D(const BasisFunc2D& bf)
	{
		// Copy standard
		this->_min1 = bf._min1;
		this->_max1 = bf._max1;
		this->_n1 = bf._n1;
		this->_delta1 = bf._delta1;
		this->_min2 = bf._min2;
		this->_max2 = bf._max2;
		this->_n2 = bf._n2;
		this->_delta2 = bf._delta2;

	    // Allocate then copy soln, nus
	    this->_vals = new double*[bf._n1];
		for (int i=0; i<bf._n1; i++) {
			this->_vals[i] = new double[bf._n2];
		    *(this->_vals[i]) = *(bf._vals[i]);
		};
	    this->_valsb = new double*[bf._n1+2];
		for (int i=0; i<bf._n1+2; i++) {
			this->_valsb[i] = new double[bf._n2+2];
		    *(this->_valsb[i]) = *(bf._valsb[i]);
		};
	};
	BasisFunc2D& BasisFunc2D::operator=(const BasisFunc2D& bf)
	{
		if (this != &bf)  
		{
			// Copy standard
			this->_min1 = bf._min1;
			this->_max1 = bf._max1;
			this->_n1 = bf._n1;
			this->_delta1 = bf._delta1;
			this->_min2 = bf._min2;
			this->_max2 = bf._max2;
			this->_n2 = bf._n2;
			this->_delta2 = bf._delta2;

			// Free the existing resource
			for (int i=0; i<this->_n1; i++) {
				safeDelArr(this->_vals[i]);  
			};
			safeDelArr(this->_vals);
			for (int i=0; i<this->_n1+2; i++) {
				safeDelArr(this->_valsb[i]);  
			};
			safeDelArr(this->_valsb);  

		    // Allocate then copy soln, nus
		    this->_vals = new double*[bf._n1];
			for (int i=0; i<bf._n1; i++) {
				this->_vals[i] = new double[bf._n2];
			    *(this->_vals[i]) = *(bf._vals[i]);
			};
		    this->_valsb = new double*[bf._n1+2];
			for (int i=0; i<bf._n1+2; i++) {
				this->_valsb[i] = new double[bf._n2+2];
			    *(this->_valsb[i]) = *(bf._valsb[i]);
			};
		};
		return *this;
	};
	BasisFunc2D::~BasisFunc2D() {
		for (int i=0; i<this->_n1; i++) {
			safeDelArr(this->_vals[i]);  
		};
		safeDelArr(this->_vals); 
		for (int i=0; i<this->_n1+2; i++) {
			safeDelArr(this->_valsb[i]);  
		};
		safeDelArr(this->_valsb);   
	};

	/********************
	Get values
	********************/

	double BasisFunc2D::get_val(double nu1, double nu2) 
	{
		return _bicubic_interp(nu1, nu2);
	};
	double BasisFunc2D::get_deriv(double nu1, double nu2, bool deriv_x, bool deriv_y) {
		return _bicubic_interp(nu1, nu2, deriv_x, deriv_y);
	};

	/********************
	Write lattice to a file
	********************/
	
	void BasisFunc2D::write_to_file(std::string fname) {
		_write_to_file(fname);
	};
	/*
	void BasisFunc2D::write_to_file(std::string fname, int npts_x, int npts_y) {
		_write_to_file(fname,npts_x,npts_y);
	};
	*/

	/********************
	IO
	********************/

	std::ostream& operator<<(std::ostream& os,  BasisFunc2D& bf)
	{
		for (int i=0; i<bf._n1; i++) 
		{ 
			for (int j=0; j<bf._n2; j++)
			{
				os << bf._min1 + i*bf._delta1 << " " << bf._min2 + j*bf._delta2 << " " << bf._vals[i][j];
				if (i != bf._n1-1 || j != bf._n2-1) {
					os << "\n";
				};
			};
		};
		return os;
	};

	/********************
	Update the values
	********************/

	void BasisFunc2D::update(double **dvals) {
		for (int i=0; i<_n1; i++) {
			for (int j=0; j<_n2; j++) {
				_vals[i][j] += dvals[i][j];
			};
		};

		// Update the bounding box
		_update_bounding_box();
	};

	/****************************************
	BasisFunc2D - PRIVATE
	****************************************/

	// For interpolation, see:
	// http://www.paulinternet.nl/?page=bicubic

	/********************
	Update the bounding box
	********************/
	void BasisFunc2D::_update_bounding_box() {
		for (int i1=0; i1<_n1+2; i1++) {
			for (int i2=0; i2<_n2+2; i2++) {
				if (i1 == 0) {
					if (i2 == 0) {
						_valsb[i1][i2] = 2.*_vals[i1][i2] - _vals[i1+1][i2+1];
					} else if (i2 == _n2+1) {
						_valsb[i1][i2] = 2.*_vals[i1][i2-2] - _vals[i1+1][i2-3];
					} else {
						_valsb[i1][i2] = 2.*_vals[i1][i2-1] - _vals[i1+1][i2-1];
					};
				} else if (i1 == _n1+1) {
					if (i2 == 0) {
						_valsb[i1][i2] = 2.*_vals[i1-2][i2] - _vals[i1-3][i2+1];
					} else if (i2 == _n2+1) {
						_valsb[i1][i2] = 2.*_vals[i1-2][i2-2] - _vals[i1-3][i2-3];
					} else {
						_valsb[i1][i2] = 2.*_vals[i1-2][i2-1] - _vals[i1-3][i2-1];
					};
				} else {
					if (i2 == 0) {
						_valsb[i1][i2] = 2.*_vals[i1-1][i2] - _vals[i1-1][i2+1];
					} else if (i2 == _n2+1) {
						_valsb[i1][i2] = 2.*_vals[i1-1][i2-2] - _vals[i1-1][i2-3];
					} else {
						_valsb[i1][i2] = _vals[i1-1][i2-1];
					};
				};
			};
		};
	};

	/********************
	Get lower indexes in bounding box
	********************/

	std::pair<int,int> BasisFunc2D::_get_idx_interval(double x, double y) {
		if (x < _min1 || x > _max1 || y < _min2 || y > _max2) {
			std::cerr << "ERROR - outside interval." << std::endl; 
			exit(EXIT_FAILURE); 
		};		
		int i = (x - _min1)/_delta1;
		if (i == _n1 - 1) { i--; };
		int j = (y - _min2)/_delta2;
		if (j == _n2 - 1) { j--; };

		// Reindex to bounding box
		i++;
		j++;

		return std::make_pair(i,j);
	};

	// Interpolate a given 4 points
	double BasisFunc2D::_cubic_interp(double *p, double x, int dim)
	{
		double delta,max;
		if (dim == 0) {
			delta = _delta1;
			max = _max1;
		} else if (dim == 1) {
			delta = _delta2;
			max = _max2;
		};
		double f = fmod(x,_delta1)/_delta1;
		if (f > 1.0-1E-6) { f=0.0; }; // always use 0 rather than 1
		if (_max1 - x < 1E-6) { f=1.0; }; // very end
		return p[1] + 0.5 * f*(p[2] - p[0] + f*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + f*(3.0*(p[1] - p[2]) + p[3] - p[0])));
	};

	// Get the deriv of 4 points
	double BasisFunc2D::_get_deriv(double *p, double x, int dim)
	{
		double delta,max;
		if (dim == 0) {
			delta = _delta1;
			max = _max1;
		} else if (dim == 1) {
			delta = _delta2;
			max = _max2;
		};
		double f = fmod(x,_delta1)/_delta1;
		if (f > 1.0-1E-6) { f=0.0; }; // always use 0 rather than 1
		if (_max1 - x < 1E-6) { f=1.0; }; // very end
		return 0.5 * (p[2] - p[0]) + 2.0 * f * (p[0] - 2.5 * p[1] + 2.0 * p[2] - 0.5 * p[3]) + 3.0 * f * f *( - 0.5 * p[0] + 1.5 * p[1] - 1.5 * p[2] + 0.5 * p[3]);
	};

	// Bicubic interpolation
	double BasisFunc2D::_bicubic_interp(double x, double y, bool deriv_x, bool deriv_y)
	{
		// Get the bounding idxs in the bounding box
		std::pair<int,int> idxs = _get_idx_interval(x,y);
		int i = idxs.first, j = idxs.second;

		// Store interp vals
		double *arr = new double[4];		

		// The vals to interpolate
		double *p = new double[4];

		// i-1
		p[0] = _valsb[i-1][j-1];
		p[1] = _valsb[i-1][j];
		p[2] = _valsb[i-1][j+1];
		p[3] = _valsb[i-1][j+2];
		if (deriv_y) {
			arr[0] = _get_deriv(p,y,1);
		} else {
			arr[0] = _cubic_interp(p,y,1);			
		};
		// i
		p[0] = _valsb[i][j-1];
		p[1] = _valsb[i][j];
		p[2] = _valsb[i][j+1];
		p[3] = _valsb[i][j+2];
		if (deriv_y) {
			arr[1] = _get_deriv(p,y,1);
		} else {
			arr[1] = _cubic_interp(p,y,1);			
		};
		// i+1
		p[0] = _valsb[i+1][j-1];
		p[1] = _valsb[i+1][j];
		p[2] = _valsb[i+1][j+1];
		p[3] = _valsb[i+1][j+2];
		if (deriv_y) {
			arr[2] = _get_deriv(p,y,1);
		} else {
			arr[2] = _cubic_interp(p,y,1);			
		};
		// i+2
		p[0] = _valsb[i+2][j-1];
		p[1] = _valsb[i+2][j];
		p[2] = _valsb[i+2][j+1];
		p[3] = _valsb[i+2][j+2];
		if (deriv_y) {
			arr[3] = _get_deriv(p,y,1);
		} else {
			arr[3] = _cubic_interp(p,y,1);			
		};

		// Interpolate one last time
		if (deriv_x) {
			return _get_deriv(arr,x,0);
		} else {
			return _cubic_interp(arr,x,0);
		};
	};

	/********************
	Write to a file
	********************/

	void BasisFunc2D::_write_to_file(std::string fname) 
	{
		std::ofstream f;
		f.open (fname);
		for (int i=0; i<_n1; i++) {
			for (int j=0; j<_n2; j++) {
				f << _min1+i*_delta1 << " " << _min2+j*_delta2 << " " << _vals[i][j];
				if (i != _n1-1 || j != _n2-1) { f << "\n";};
			};
		};
		f.close();
	};

	/*
	void BasisFunc2D::_write_to_file(std::string fname, int npts_x, int npts_y)
	{
		std::ofstream f;
		f.open (fname);
		double d1 = (_max1-_min1) / (npts_x-1);
		double d2 = (_max2-_min2) / (npts_y-1);
		for (int i=0; i<npts_x; i++) {
			for (int j=0; j<npts_y; j++) {
				f << _min1+i*d1 << " " << _min2+j*d2 << " " << _cubic_interp(_min1+i*d1, _min2+j*d2);
				if (i != npts_x-1 || j != npts_y-1) { f << "\n";};
			};
		};
		f.close();
	};
	*/

};