#ifndef IOSTREAM_h
#define IOSTREAM_h
#include <iostream>
#endif

#ifndef VECTOR_h
#define VECTOR_h
#include <vector>
#endif

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	BasisFunc
	****************************************/

	class BasisFunc
	{
	private:

		// Derivative
		double* _deriv;

		// Vals
		double* _vals;

		// Get the two points left and right of a given point
		double *_get_interval(double *s, double x);

		// Interpolation functions
		double _cubic_interp(double *s, double x);

		// Deriv of soln at pt
		double _get_deriv(double x);

		// Write function
		void _write_to_file(double *s, std::string fname);
		void _write_to_file(double *s, std::string fname, int npts);

		// Range
		double _min;
		double _max;

		// Number increments
		int _n;

		// Increment size
		double _delta;

		// Calculate derivative at all points and store
		void _update_derivs();

	public:

		// Constructor
		BasisFunc(double nu_min, double nu_max, int n_nu);
		BasisFunc(const BasisFunc& bf);
		BasisFunc& operator=(const BasisFunc& bf);
		~BasisFunc();

		// Get values, if they are in the lattice
		double get_val(double nu);
		double get_deriv(double nu);

		// Write deriv to a file
		void write_to_file(std::string fname);
		void write_to_file(std::string fname, int npts);
		void write_deriv_to_file(std::string fname);
		void write_deriv_to_file(std::string fname, int npts);

		// To write private vars, declare as friend
		friend std::ostream& operator<<(std::ostream& os, BasisFunc& bf);

		// Accessor
		double& operator[](std::size_t idx);

		// Update the values
		void update(double *dvals);

		// Update the values, smoothing delta F before
		void update_smoothed(double *dvals);

		// Test: fill with a sin func
		void test_fill_sin();
	};

};