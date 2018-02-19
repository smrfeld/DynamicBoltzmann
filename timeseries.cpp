
#include <iostream>

#ifndef TIMESERIES_h
#define TIMESERIES_h
#include "timeseries.hpp"
#endif

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Timeseries
	****************************************/

	/********************
	Constructors
	********************/

	Timeseries::Timeseries(double tmin, double tmax, int n_t)
	{
		this->_tmin = tmin;
		this->_tmax = tmax;
		this->_n_t = n_t;
		this->_dt = (tmax-tmin) / n_t;
		this->_soln = new double[n_t];
	};
	Timeseries::Timeseries(const Timeseries& t)
	{
		// Copy standard
		this->_tmin = t._tmin;
		this->_tmax = t._tmax;
		this->_n_t = t._n_t;
		this->_dt = t._dt;

	    // Allocate memory for the new soln
	    this->_soln = new double[t._n_t];
	    // Copy
	    *(this->_soln) = *(t._soln);
	};
	Timeseries& Timeseries::operator=(const Timeseries& t) 
	{
		if (this != &t)  
		{  
			// Copy standard
			this->_tmin = t._tmin;
			this->_tmax = t._tmax;
			this->_n_t = t._n_t;
			this->_dt = t._dt;

			 // Free the existing resource.  
			 delete this->_soln;  
				
			// Allocate memory for the new soln
			this->_soln = new double[t._n_t];
			// Copy
			*(this->_soln) = *(t._soln);

		};
		return *this;
	};
	Timeseries::~Timeseries() {
		delete this->_soln;
	};


};