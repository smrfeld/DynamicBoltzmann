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
	BasisFunc1D
	****************************************/

	class BasisFunc1D
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
		BasisFunc1D(double nu_min, double nu_max, int n_nu);
		BasisFunc1D(const BasisFunc1D& bf);
		BasisFunc1D& operator=(const BasisFunc1D& bf);
		~BasisFunc1D();

		// Get values, if they are in the lattice
		double get_val(double nu);
		double get_deriv(double nu);

		// Write deriv to a file
		void write_to_file(std::string fname);
		void write_to_file(std::string fname, int npts);
		void write_deriv_to_file(std::string fname);
		void write_deriv_to_file(std::string fname, int npts);

		// To write private vars, declare as friend
		friend std::ostream& operator<<(std::ostream& os, BasisFunc1D& bf);

		// Accessor
		double& operator[](std::size_t idx);

		// Update the values
		void update(double *dvals);

		// Update the values, smoothing delta F before
		void update_smoothed(double *dvals);

		// Test: fill with a sin func
		void test_fill_sin();
	};


	/****************************************
	BasisFunc2D
	****************************************/

	class BasisFunc2D
	{
	private:

		// Vals
		double** _vals;

		// Bounded vals
		double** _valsb;

		// Get the two points left and right of a given point
		double *_get_interval(double **s, double x, int dim);

		// Interpolation functions
		double _cubic_interp(double *p, double x, int dim);
		double _bicubic_interp(double x, double y, bool deriv_x=false, bool deriv_y=false);

		// Deriv of soln at pt
		double _get_deriv(double *p, double x, int dim);

		// Write function
		void _write_to_file(std::string fname);
		// void _write_to_file(std::string fname, int npts_x, int npts_y);

		// Range
		double _min1,_max1,_min2,_max2;

		// Number increments
		int _n1,_n2;

		// Increment size
		double _delta1,_delta2;

		// Get lower indexes in bounding box
		std::pair<int,int> _get_idx_interval(double x, double y);

		/********************
		Update the bounding box
		********************/
		void _update_bounding_box();

	public:

		// Constructor
		BasisFunc2D(double nu_min1, double nu_max1, int n_nu1, double nu_min2, double nu_max2, int n_nu2);
		BasisFunc2D(const BasisFunc2D& bf);
		BasisFunc2D& operator=(const BasisFunc2D& bf);
		~BasisFunc2D();

		// Get values, if they are in the lattice
		double get_val(double nu1, double nu2);
		double get_deriv(double nu1, double nu2, bool deriv_x, bool deriv_y);

		// Write deriv to a file
		void write_to_file(std::string fname);
		// void write_to_file(std::string fname, int npts_x, int npts_y);

		// To write private vars, declare as friend
		friend std::ostream& operator<<(std::ostream& os, BasisFunc2D& bf);

		// Update the values
		void update(double **dvals);
	};

};