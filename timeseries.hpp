#ifndef IOSTREAM_h
#define IOSTREAM_h
#include <iostream>
#endif

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Timeseries
	****************************************/

	class Timeseries
	{
	private:

		// Range in time
		double _t_min;
		double _t_max;

		// Number timesteps
		int _n_t;

		// Timestep
		double _dt;

		// Grid to store solution
		// Pointer so size can be dynamically allocated by constructor
		double* _tvals;
		double* _soln;

	public:

		// Constructor
		Timeseries(double t_min, double t_max, int n_t);
		Timeseries(const Timeseries& t);
		Timeseries& operator=(const Timeseries& t);
		~Timeseries();

		// To write private vars, declare as friend
		friend std::ostream& operator<<(std::ostream& os, Timeseries& t);
	};
};