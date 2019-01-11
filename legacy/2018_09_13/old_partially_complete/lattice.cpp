#include "lattice.hpp"
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
	Lattice
	****************************************/

	// Constructor
	Lattice::Lattice(double min, double max, int n)
	{
		this->min = min;
		this->max = max;
		this->n = n;
		this->delta = (max - min) / (n - 1);

		// Zero by default
	    this->vals = new double[n];
		std::fill_n(this->vals, n, 0.0);
	};
	Lattice::Lattice(const Lattice& l)
	{
		// Copy standard
		this->min = l.min;
		this->max = l.max;
		this->n = l.n;
		this->delta = l.delta;

	    // Allocate then copy soln, nus
	    this->vals = new double[l.n];
	    *(this->vals) = *(l.vals);
	};
	Lattice& Lattice::operator=(const Lattice& l) 
	{
		if (this != &l)  
		{  
			// Copy standard
			this->min = l.min;
			this->max = l.max;
			this->n = l.n;
			this->delta = l.delta;

			// Free the existing resource.  
			safeDelArr(this->vals);  

		    // Allocate then copy soln, nus
		    this->vals = new double[l.n];
		    *(this->vals) = *(l.vals);
		};
		return *this;
	};
	Lattice::~Lattice() {
		safeDelArr(this->vals);  
	};

	/********************
	IO
	********************/

	std::ostream& operator<<(std::ostream& os,  Lattice& l)
	{
		for (int i=0; i<l.n; i++) 
		{ 
			os << l.min + i*l.delta << " " << l[i];
			if (i != l.n-1) {
				os << "\n";
			};
		};
		return os;
	};

	/********************
	Write lattice to a file
	********************/

	void Lattice::write_to_file(std::string fname) 
	{
		_write_to_file(vals,fname);
	};

	/********************
	Accessor
	********************/

	// Accessor
	double& Lattice::operator[](std::size_t idx) { return vals[idx]; };
	
	/****************************************
	Lattice - PROTECTED
	****************************************/

	/********************
	Write to file
	********************/

	void Lattice::_write_to_file(double *s, std::string fname) 
	{
		std::ofstream f;
		f.open (fname);
		for (int i=0; i<n; i++) {
			f << min+i*delta << " " << s[i];
			if (i != n) { f << "\n";};
		};
		f.close();
	};

	/****************************************
	BasisFunc
	****************************************/

	// Constructor
	BasisFunc::BasisFunc(double nu_min, double nu_max, int n_nu) : Lattice(nu_min,nu_max,n_nu) 
	{
		// Check min no pts is satisfied
		if (n_nu <= 4) { 
			std::cerr << "ERROR - n_nu > 4 required." << std::endl; 
			exit(EXIT_FAILURE); 
		};

		this->_deriv = new double[n_nu];
		// All zero by default
		std::fill_n(this->_deriv, n_nu, 0.0);
	};
	BasisFunc::BasisFunc(const BasisFunc& bf) : Lattice(bf)
	{
	    // Allocate then copy soln, nus
	    this->_deriv = new double[bf.n];
	    *(this->_deriv) = *(bf._deriv);
	};
	BasisFunc& BasisFunc::operator=(const BasisFunc& bf)
	{
		if (this != &bf)  
		{  
	        Lattice::operator=(bf);

			// Free the existing resource.  
			safeDelArr(this->_deriv);  

		    // Allocate then copy soln, nus
		    this->_deriv = new double[bf.n];
		    *(this->_deriv) = *(bf._deriv);
		};
		return *this;
	};
	BasisFunc::~BasisFunc() {
		safeDelArr(this->_deriv);  
	};

	/********************
	Calculate derivative at all points and store
	********************/

	void BasisFunc::update_derivs() {
		for (int i=0; i<n; i++) {
			_deriv[i] = _get_deriv(min+i*delta);
		};
	};

	/********************
	Get values
	********************/

	double BasisFunc::get_val(double nu) 
	{
		return _cubic_interp(this->vals, nu);
	};
	double BasisFunc::get_deriv(double nu) {
		return _cubic_interp(this->_deriv, nu);
	};

	/********************
	Write lattice to a file
	********************/
	
	void BasisFunc::write_to_file(std::string fname) {
		Lattice::_write_to_file(vals,fname);
	};
	void BasisFunc::write_to_file(std::string fname, int npts) {
		_write_to_file(vals,fname,npts);
	};
	void BasisFunc::write_deriv_to_file(std::string fname) {
		Lattice::_write_to_file(_deriv,fname);
	};
	void BasisFunc::write_deriv_to_file(std::string fname, int npts) {
		_write_to_file(_deriv,fname,npts);
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
		if (x < min || x > max) { 
			std::cerr << "ERROR - outside interval." << std::endl; 
			exit(EXIT_FAILURE); 
		};		
		int i = (x - min)/delta;
		if (i == n - 1) { i--; };
		static double p[4];
		p[1] = s[i];
		p[2] = s[i+1];

		if (i == 0) {
			// Between the first two points - make line for first point
			p[0] = 2.*s[i] - s[i+1];
			p[3] = s[i+2];
		} else if (i == n - 2) {
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
		double f = fmod(x,delta)/delta;
		if (f > 1.0-1E-6) { f=0.0; }; // always use 0 rather than 1
		if (max - x < 1E-6) { f=1.0; }; // very end
		return p[1] + 0.5 * f*(p[2] - p[0] + f*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + f*(3.0*(p[1] - p[2]) + p[3] - p[0])));
	};

	double BasisFunc::_get_deriv(double x)
	{
		double *p = _get_interval(vals,x);
		double f = fmod(x,delta)/delta;
		if (f > 1.0-1E-6) { f=0.0; }; // always use 0 rather than 1
		if (max - x < 1E-6) { f=1.0; }; // very end
		return 0.5 * (p[2] - p[0]) + 2.0 * f * (p[0] - 2.5 * p[1] + 2.0 * p[2] - 0.5 * p[3]) + 3.0 * f * f *( - 0.5 * p[0] + 1.5 * p[1] - 1.5 * p[2] + 0.5 * p[3]);
	};

	/********************
	Write lattice to a file
	********************/

	void BasisFunc::_write_to_file(double *s, std::string fname, int npts)
	{
		std::ofstream f;
		f.open (fname);
		double d = (max-min) / (npts-1);
		for (int i=0; i<npts; i++) {
			f << min + i*d << " " << _cubic_interp(s, min + i*d);
			if (i != npts) { f << "\n";};
		};
		f.close();
	};


};