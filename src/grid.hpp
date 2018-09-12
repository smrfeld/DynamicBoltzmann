#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

#ifndef STRING_H
#define STRING_H
#include <string>
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	Grid
	****************************************/

	class Grid {

	private:

		// Name
		std::string _name;

		// Size
		int _n;
		// Min/max/increment
		double _min;
		double _max;
		double _delta;
		// Grid of vals
		double *_grid;

		// Copy/clean up
		void _copy(const Grid& other);
		void _clean_up();

	public:

		// Constructors
		Grid(std::string name, double min, double max, int n);
		Grid(const Grid& other);
		Grid & operator=(const Grid& other);
		~Grid();

		// Getters/setters
		std::string name() const;
		double delta() const;
		int n() const;
		double get_by_idx(int i) const;

		// Set n
		void set_n(int n);

		// Get indexes surrounding a point
		// ie point is between i and i+1 where i is returned
		int surrounding_idxs(double x) const; 

		// Get fraction of a point between successive points
		double frac_between(double x) const;
		// Second optional specification: the return of the surrounding idxs
		double frac_between(double x, int i) const;

		// Check if a given point is in the grid
		bool in_grid(double x) const;

		// Print grid range
		void print_grid_range() const;

		// Write the grid into an ofstream
		void write_grid(std::string fname) const;

		// Test: values if a sin function were defined on the grid
		std::vector<double> test_sin() const;
		std::vector<double> test_cos() const;
	};

};