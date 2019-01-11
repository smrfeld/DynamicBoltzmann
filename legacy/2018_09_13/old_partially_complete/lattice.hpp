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
	Lattice2D
	****************************************/

	class Lattice2D
	{
	protected:

		// Write functions
		void _write_to_file(double *s, std::string fname);

	public:

		// Ranges
		double min1,max1,min2,max2;

		// Number increments
		int n1,n2;

		// Increment size
		double delta1,delta2;

		// Values
		// Pointer so size can be dynamically allocated by constructor
		double **vals;

		// Constructor
		Lattice2D(double min1, double max1, int n1, double min2, double max2, int n2);
		Lattice2D(const Lattice2D& l);
		Lattice2D& operator=(const Lattice2D& l);
		~Lattice2D();

		// To write private vars, declare as friend
		friend std::ostream& operator<<(std::ostream& os, Lattice2D& l);

		// Write lattice to a file
		void write_to_file(std::string fname);

		// Accessor
		double*& operator[](std::size_t idx);
	};

	/****************************************
	Lattice
	****************************************/

	class Lattice
	{
	protected:

		// Write functions
		void _write_to_file(double *s, std::string fname);

	public:

		// Range
		double min;
		double max;

		// Number increments
		int n;

		// Increment size
		double delta;

		// Values
		// Pointer so size can be dynamically allocated by constructor
		double *vals;

		// Constructor
		Lattice(double min, double max, int n);
		Lattice(const Lattice& l);
		Lattice& operator=(const Lattice& l);
		~Lattice();

		// To write private vars, declare as friend
		friend std::ostream& operator<<(std::ostream& os, Lattice& l);

		// Write lattice to a file
		void write_to_file(std::string fname);

		// Accessor
		double& operator[](std::size_t idx);
	};

	/****************************************
	BasisFunc
	****************************************/

	class BasisFunc : public Lattice
	{
	private:

		// Derivative
		double* _deriv;

		// Get the two points left and right of a given point
		double *_get_interval(double *s, double x);

		// Interpolation functions
		double _cubic_interp(double *s, double x);

		// Deriv of soln at pt
		double _get_deriv(double x);

		// Write function
		void _write_to_file(double *s, std::string fname, int npts);

	public:

		// Constructor
		BasisFunc(double nu_min, double nu_max, int n_nu);
		BasisFunc(const BasisFunc& bf);
		BasisFunc& operator=(const BasisFunc& bf);
		~BasisFunc();

		// Calculate derivative at all points and store
		void update_derivs();

		// Get values, if they are in the lattice
		double get_val(double nu);
		double get_deriv(double nu);

		// Write deriv to a file
		void write_to_file(std::string fname);
		void write_to_file(std::string fname, int npts);
		void write_deriv_to_file(std::string fname);
		void write_deriv_to_file(std::string fname, int npts);

	};

};