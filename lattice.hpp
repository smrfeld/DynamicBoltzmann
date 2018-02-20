#ifndef IOSTREAM_h
#define IOSTREAM_h
#include <iostream>
#endif

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Lattice
	****************************************/

	class Lattice
	{
	private:

		// Range in 1d
		double _nu_min;
		double _nu_max;
		
		// Number increments
		int _n_nu;

		// Increment size
		double _dnu;

		// Grid to store solution
		// Pointer so size can be dynamically allocated by constructor
		double* _nuvals;
		double* _soln;
		double* _soln_deriv;

		// Get the two points left and right of a given point
		double *_get_interval(double *s, double x);

		// Interpolation functions
		double _cubic_interp(double *s, double x);

		// Deriv of soln at pt
		double _deriv(double x);

		// Calculate derivative at all points and store
		void _diff_soln();

		// Write functions
		void _write_to_file(double *s, std::string fname);
		void _write_to_file(double *s, std::string fname, int npts);

	public:

		// Constructor
		Lattice(double nu_min, double nu_max, int n_nu);
		Lattice(const Lattice& l);
		Lattice& operator=(const Lattice& l);
		~Lattice();

		// To write private vars, declare as friend
		friend std::ostream& operator<<(std::ostream& os, Lattice& l);

		// tmp
		void test_fill();

		// Write lattice to a file
		void write_to_file(std::string fname);
		void write_to_file(std::string fname, int npts);
		void write_deriv_to_file(std::string fname);
		void write_deriv_to_file(std::string fname, int npts);

	};

};