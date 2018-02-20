#include <iostream>
#include <algorithm>
#include "timeseries.hpp"
#include "general.hpp"

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

	Timeseries::Timeseries(double t_min, double t_max, int n_t)
	{
		this->_t_min = t_min;
		this->_t_max = t_max;
		this->_n_t = n_t;
		this->_dt = (t_max-t_min) / (n_t-1);

		this->_soln = new double[n_t];
		// All zero by default
		std::fill_n(this->_soln, n_t, 0.0);

		// Time values
		this->_tvals = new double[n_t];
		for (int i=0; i<n_t; i++) {
			this->_tvals[i] = t_min + i*this->_dt;
		};
	};
	Timeseries::Timeseries(const Timeseries& t)
	{
		// Copy standard
		this->_t_min = t._t_min;
		this->_t_max = t._t_max;
		this->_n_t = t._n_t;
		this->_dt = t._dt;

	    // Allocate then copy soln, times
	    this->_soln = new double[t._n_t];
	    *(this->_soln) = *(t._soln);
	    this->_tvals = new double[t._n_t];
	    *(this->_tvals) = *(t._tvals);
	};
	Timeseries& Timeseries::operator=(const Timeseries& t) 
	{
		if (this != &t)  
		{  
			// Copy standard
			this->_t_min = t._t_min;
			this->_t_max = t._t_max;
			this->_n_t = t._n_t;
			this->_dt = t._dt;

			// Free the existing resource.  
			safeDelArr(this->_soln);  
			safeDelArr(this->_tvals);

		    // Allocate then copy soln, times
			this->_soln = new double[t._n_t];
			*(this->_soln) = *(t._soln);
		    this->_tvals = new double[t._n_t];
		    *(this->_tvals) = *(t._tvals);
		};
		return *this;
	};
	Timeseries::~Timeseries() {
		safeDelArr(this->_soln);  
		safeDelArr(this->_tvals);
	};

	/********************
	IO
	********************/

	std::ostream& operator<<(std::ostream& os,  Timeseries& t)
	{
		for (int i=0; i<t._n_t; i++) 
		{ 
			os << t._tvals[i] << " " << t._soln[i];
			if (i != t._n_t-1) {
				os << "\n";
			};
		};
		return os;
	};

	/****************************************
	Timeseries - PRIVATE
	****************************************/




};