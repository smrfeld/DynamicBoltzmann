#include <string>
#include <map>

#ifndef FWDS_IXN_PARAM_H
#define FWDS_IXN_PARAM_H
#include "fwds/fwds_ixn_param.hpp"
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

		// Name
		std::string _name;

		// The ixn param
		Iptr _ixn_param;

		// Values
		// Timepoints = timesteps + 1
		int _no_timepoints;
		int _no_timesteps;
		double *_vals;
		int _zero_end_cond_timepoint;

		// Internal copy func/clean up
		void _clean_up();
		void _copy(const Adjoint& other);
		void _move(Adjoint &other);

	public:

		/********************
		Constructor
		********************/

		Adjoint(std::string name, Iptr ixn_param);
		Adjoint(const Adjoint& other);
		Adjoint& operator=(const Adjoint& other);
		Adjoint(Adjoint&& other);
		Adjoint& operator=(Adjoint&& other);
		~Adjoint();

		/********************
		Check setup
		********************/

		void check_setup() const;

		/********************
		Timesteps
		********************/

		int get_no_timesteps() const;
		void set_no_timesteps(int no_timesteps);

		/********************
		End cond
		********************/

		int get_zero_end_cond_timepoint() const;
		void set_zero_end_cond_timepoint(int timepoint);

		/********************
		Getters (general)
		********************/

		std::string get_name() const;
		Iptr get_ixn_param() const;

		/********************
		Value
		********************/

		double get_val_at_timepoint(int timepoint) const;

		/********************
		Solve diff eq
		********************/

		void solve_diff_eq_at_timepoint_to_minus_one(int timepoint, double dt, bool l2_mode=false, const std::map<Iptr,double> &l2_lambda = std::map<Iptr,double>(), const std::map<Iptr,double> &l2_center = std::map<Iptr,double>());

		/********************
		Reset to zero
		********************/

		void reset_to_zero();

		/********************
		Write to file
		********************/

		void write_to_file(std::string fname) const;
	};

};