#include <string>
#include <map>
#include <vector>

#include "fwds/fwds_ixn_param_traj.hpp"

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
		ITptr _ixn_param_traj;

		// Values
		// Timepoints = timesteps + 1
		int _no_timepoints;
		int _no_timesteps;
        std::vector<double> _vals;
		int _timepoint_zero_end_cond;

		// Internal copy func/clean up
		void _clean_up();
		void _copy(const Adjoint& other);
		void _move(Adjoint &other);

	public:

        // ***************
        // MARK: - Constructor
        // ***************

		Adjoint(std::string name, ITptr ixn_param_traj);
		Adjoint(const Adjoint& other);
		Adjoint& operator=(const Adjoint& other);
		Adjoint(Adjoint&& other);
		Adjoint& operator=(Adjoint&& other);
		~Adjoint();

        // ***************
        // MARK: - Timesteps
        // ***************
        
		int get_no_timesteps() const;
		void set_no_timesteps(int no_timesteps);

        // ***************
        // MARK: - End condition
        // ***************
        
		int get_timepoint_zero_end_cond() const;
		void set_timepoint_zero_end_cond(int timepoint);

        // ***************
        // MARK: - Getters (general)
        // ***************

		std::string get_name() const;
		ITptr get_ixn_param_traj() const;

        // ***************
        // MARK: - Get value at timepoint
        // ***************

		double get_val_at_timepoint(int timepoint) const;

        // ***************
        // MARK: - Solve diff eq
        // ***************

		void solve_diff_eq_at_timepoint_to_minus_one(int timepoint, double dt, bool l2_mode=false, const std::map<ITptr,double> &l2_lambda = std::map<ITptr,double>(), const std::map<ITptr,double> &l2_center = std::map<ITptr,double>());

        // ***************
        // MARK: - Write to file
        // ***************
        
		void write_to_file(std::string fname) const;
	};

};
