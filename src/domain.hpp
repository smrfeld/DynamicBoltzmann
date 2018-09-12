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
	Domain1D
	****************************************/

	class Domain1D {

	private:

		// Name
		std::string _name;

		// Size
		int _no_pts;

		// Min/max/increment
		double _min;
		double _max;
		double _delta;

		// Copy/clean up
		void _copy(const Domain1D& other);
		void _clean_up();

	public:

		/********************
		Constructor
		********************/

		Domain1D(std::string name, double min, double max, int no_pts);
		Domain1D(const Domain1D& other);
		Domain1D(Domain1D&& other);
		Domain1D& operator=(const Domain1D& other);
		Domain1D& operator=(Domain1D&& other);
		~Domain1D();

		/********************
		Getters
		********************/

		std::string get_name() const;
		double get_delta() const;
		int get_no_pts() const;

		/********************
		Value of a point
		********************/

		double get_by_idx(int i) const;

		/********************
		Check if point is in domain
		********************/

		bool check_in_domain(double x) const;

		/********************
		Resize
		********************/

		void resize_no_pts_at_fixed_endpoints(int no_pts);
		void resize_no_pts_at_fixed_spacing(int no_pts);

		/********************
		Get idxs/fractions for points
		********************/

		// Get indexes surrounding a point
		// ie point is between i and i+1 where i is returned
		int get_surrounding_idxs(double x) const; 

		// Get fraction of a point between successive points
		double get_frac_between(double x) const;
		// Second optional specification: the return of the surrounding idxs
		double get_frac_between(double x, int i) const;

		/********************
		I/O
		********************/

		// Print domain range
		void print_domain_range() const;

		// Write the domain into an ofstream
		void write_domain(std::string fname) const;
	};

};