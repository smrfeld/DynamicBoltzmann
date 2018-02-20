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
	Lattice::Lattice(double nu_min, double nu_max, int n_nu)
	{
		// Check min no pts is satisfied
		if (n_nu <= 4) { 
			std::cerr << "ERROR - n_nu > 4 required." << std::endl; 
			exit(EXIT_FAILURE); 
		};

		this->_nu_min = nu_min;
		this->_nu_max = nu_max;
		this->_n_nu = n_nu;
		this->_dnu = (nu_max - nu_min) / (n_nu - 1);

		this->_soln = new double[n_nu];
		// All zero by default
		std::fill_n(this->_soln, n_nu, 0.0);

		this->_soln_deriv = new double[n_nu];
		// All zero by default
		std::fill_n(this->_soln_deriv, n_nu, 0.0);

		// Time values
		this->_nuvals = new double[n_nu];
		for (int i=0; i<n_nu; i++) {
			this->_nuvals[i] = nu_min + i*this->_dnu;
		};
	};
	Lattice::Lattice(const Lattice& l)
	{
		// Copy standard
		this->_nu_min = l._nu_min;
		this->_nu_max = l._nu_max;
		this->_n_nu = l._n_nu;
		this->_dnu = l._dnu;

	    // Allocate then copy soln, nus
	    this->_soln = new double[l._n_nu];
	    *(this->_soln) = *(l._soln);
	    this->_soln_deriv = new double[l._n_nu];
	    *(this->_soln_deriv) = *(l._soln_deriv);
	    this->_nuvals = new double[l._n_nu];
	    *(this->_nuvals) = *(l._nuvals);
	};
	Lattice& Lattice::operator=(const Lattice& l) 
	{
		if (this != &l)  
		{  
			// Copy standard
			this->_nu_min = l._nu_min;
			this->_nu_max = l._nu_max;
			this->_n_nu = l._n_nu;
			this->_dnu = l._dnu;

			// Free the existing resource.  
			safeDelArr(this->_soln);  
			safeDelArr(this->_soln_deriv);  
			safeDelArr(this->_nuvals);

		    // Allocate then copy soln, nus
		    this->_soln = new double[l._n_nu];
		    *(this->_soln) = *(l._soln);
		    this->_soln_deriv = new double[l._n_nu];
		    *(this->_soln_deriv) = *(l._soln_deriv);
		    this->_nuvals = new double[l._n_nu];
		    *(this->_nuvals) = *(l._nuvals);
		};
		return *this;
	};
	Lattice::~Lattice() {
		safeDelArr(this->_soln);  
		safeDelArr(this->_soln_deriv);  
		safeDelArr(this->_nuvals);
	};

	/********************
	IO
	********************/

	std::ostream& operator<<(std::ostream& os,  Lattice& l)
	{
		for (int i=0; i<l._n_nu; i++) 
		{ 
			os << l._nuvals[i] << " " << l._soln[i];
			if (i != l._n_nu-1) {
				os << "\n";
			};
		};
		return os;
	};

	void Lattice::test_fill() {
		// Fill with sin
		for (int i=0; i<_n_nu; i++) {
			_soln[i] = sin(2.0*M_PI*i/(_n_nu-1))+0.1;
		};
		// Calculate deriv
		_diff_soln();
	};

	/********************
	Write lattice to a file
	********************/

	void Lattice::write_to_file(std::string fname) 
	{
		_write_to_file(_soln,fname);
	};
	void Lattice::write_to_file(std::string fname, int npts) {
		_write_to_file(_soln,fname,npts);
	};
	void Lattice::write_deriv_to_file(std::string fname) 
	{
		_write_to_file(_soln_deriv,fname);
	};
	void Lattice::write_deriv_to_file(std::string fname, int npts) {
		_write_to_file(_soln_deriv,fname,npts);
	};

	/****************************************
	Lattice - PRIVATE
	****************************************/

	// For interpolation, see:
	// http://www.paulinternet.nl/?page=bicubic

	/********************
	Get the two points left and right of a given point
	********************/

	double *Lattice::_get_interval(double *s, double x) {
		if (x < this->_nu_min || x > this->_nu_max) { 
			std::cerr << "ERROR - outside interval." << std::endl; 
			exit(EXIT_FAILURE); 
		};		
		int i = (x - this->_nu_min)/this->_dnu;
		if (i == this->_n_nu - 1) { i--; };
		static double p[4];
		p[1] = s[i];
		p[2] = s[i+1];

		if (i == 0) {
			// Between the first two points - make line for first point
			p[0] = 2.*s[i] - s[i+1];
			p[3] = s[i+2];
		} else if (i == this->_n_nu - 2) {
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
	double Lattice::_cubic_interp(double *s, double x)
	{
		double *p = _get_interval(s, x);
		double f = fmod(x,_dnu)/_dnu;
		if (f > 1.0-1E-6) { f=0.0; }; // always use 0 rather than 1
		if (_nu_max - x < 1E-6) { f=1.0; }; // very end
		return p[1] + 0.5 * f*(p[2] - p[0] + f*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + f*(3.0*(p[1] - p[2]) + p[3] - p[0])));
	};

	double Lattice::_deriv(double x)
	{
		double *p = _get_interval(_soln,x);
		double f = fmod(x,_dnu)/_dnu;
		if (f > 1.0-1E-6) { f=0.0; }; // always use 0 rather than 1
		if (_nu_max - x < 1E-6) { f=1.0; }; // very end
		return 0.5 * (p[2] - p[0]) + 2.0 * f * (p[0] - 2.5 * p[1] + 2.0 * p[2] - 0.5 * p[3]) + 3.0 * f * f *( - 0.5 * p[0] + 1.5 * p[1] - 1.5 * p[2] + 0.5 * p[3]);
	};

	/********************
	Calculate derivative at all points and store
	********************/

	void Lattice::_diff_soln() {
		for (int i=0; i<_n_nu; i++) {
			_soln_deriv[i] = _deriv(_nuvals[i]);
		};
	};

	/********************
	Write to file
	********************/

	void Lattice::_write_to_file(double *s, std::string fname) 
	{
		std::ofstream f;
		f.open (fname);
		for (int i=0; i<_n_nu; i++) {
			f << _nuvals[i] << " " << s[i];
			if (i != _n_nu) { f << "\n";};
		};
		f.close();
	};
	void Lattice::_write_to_file(double *s, std::string fname, int npts)
	{
		std::ofstream f;
		f.open (fname);
		double dnu = (_nu_max-_nu_min) / (npts-1);
		for (int i=0; i<npts; i++) {
			f << _nu_min + i*dnu << " " << _cubic_interp(s, _nu_min + i*dnu);
			if (i != npts) { f << "\n";};
		};
		f.close();
	};

};