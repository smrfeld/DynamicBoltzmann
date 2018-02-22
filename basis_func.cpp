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
	BasisFunc
	****************************************/

	// Constructor
	BasisFunc::BasisFunc(double nu_min, double nu_max, int n_nu)
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
	BasisFunc::BasisFunc(const BasisFunc& bf)
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
	BasisFunc& BasisFunc::operator=(const BasisFunc& bf)
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
	BasisFunc::~BasisFunc() {
		safeDelArr(this->_vals);
		safeDelArr(this->_deriv);  
	};

	/********************
	Get values
	********************/

	double BasisFunc::get_val(double nu) 
	{
		return _cubic_interp(_vals, nu);
	};
	double BasisFunc::get_deriv(double nu) {
		return _cubic_interp(_deriv, nu);
	};

	/********************
	Write lattice to a file
	********************/
	
	void BasisFunc::write_to_file(std::string fname) {
		_write_to_file(_vals,fname);
	};
	void BasisFunc::write_to_file(std::string fname, int npts) {
		_write_to_file(_vals,fname,npts);
	};
	void BasisFunc::write_deriv_to_file(std::string fname) {
		_write_to_file(_deriv,fname);
	};
	void BasisFunc::write_deriv_to_file(std::string fname, int npts) {
		_write_to_file(_deriv,fname,npts);
	};

	/********************
	IO
	********************/

	std::ostream& operator<<(std::ostream& os,  BasisFunc& bf)
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

	double& BasisFunc::operator[](std::size_t idx) { return _vals[idx]; };

	/********************
	Update the values
	********************/

	void BasisFunc::update(double *dvals) {
		for (int i=0; i<_n; i++) {
			_vals[i] += dvals[i];
		};

		// Update the deriv
		_update_derivs();
	};

	/********************
	Update the values, smoothing delta F before
	********************/
	void BasisFunc::update_smoothed(double *dvals)
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

	void BasisFunc::test_fill_sin() {
	    for (int i=0; i<_n; i++) {
			this->_vals[i] = sin(2.0*M_PI*i/(_n-1))+0.1;
		};
		_update_derivs();
	};

	/****************************************
	BasisFunc - PRIVATE
	****************************************/

	// For interpolation, see:
	// http://www.paulinternet.nl/?page=bicubic

	/********************
	Get the two points left and right of a given point
	********************/

	double *BasisFunc::_get_interval(double *s, double x) {
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
	double BasisFunc::_cubic_interp(double *s, double x)
	{
		double *p = _get_interval(s, x);
		double f = fmod(x,_delta)/_delta;
		if (f > 1.0-1E-6) { f=0.0; }; // always use 0 rather than 1
		if (_max - x < 1E-6) { f=1.0; }; // very end
		return p[1] + 0.5 * f*(p[2] - p[0] + f*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + f*(3.0*(p[1] - p[2]) + p[3] - p[0])));
	};

	double BasisFunc::_get_deriv(double x)
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

	void BasisFunc::_write_to_file(double *s, std::string fname) 
	{
		std::ofstream f;
		f.open (fname);
		for (int i=0; i<_n; i++) {
			f << _min+i*_delta << " " << s[i];
			if (i != _n-1) { f << "\n";};
		};
		f.close();
	};

	void BasisFunc::_write_to_file(double *s, std::string fname, int npts)
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

	void BasisFunc::_update_derivs() {
		for (int i=0; i<_n; i++) {
			_deriv[i] = _get_deriv(_min+i*_delta);
		};
	};

};