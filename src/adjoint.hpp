#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

#ifndef STRING_H
#define STRING_H
#include <string>
#endif

#ifndef MAP_H
#define MAP_H
#include <map>
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {	

	/****************************************
	Adjoint
	****************************************/

	class Adjoint {

	private:

		// The real ixn param
		IxnParamTraj *_ixn_param;

		// Name
		std::string _name;

		// No timepoints in the current traj
		int _no_timepoints_working;

		// Values
		double *_vals;

		// Internal copy func/clean up
		void _copy(const Adjoint& other);
		void _clean_up();

	public:

		// Constructor
		Adjoint(IxnParamTraj *ixn_param);
		Adjoint(const Adjoint& other);
		Adjoint& operator=(const Adjoint& other);
		Adjoint(Adjoint&& other);
		Adjoint& operator=(Adjoint&& other);
		~Adjoint();

		/********************
		Get ixn param
		********************/

		IxnParamTraj* get_ixn_param() const;

		/********************
		Solve diff eq
		********************/

		void solve_diff_eq_from_timepoint_to_timepoint_minus_one(int timepoint, double dt);
	};

};