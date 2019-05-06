#include <string>
#include <map>
#include <vector>

#include "fwds/fwds_ixn_param_traj.hpp"

/************************************
* Namespace for dblz
************************************/

namespace dblz {	

    // ***************
    // MARK: - Abstract base for the adjoint class
    // ***************
    
	class Adjoint {

	protected:

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
		virtual ~Adjoint();

        // ***************
        // MARK: - Timesteps
        // ***************
        
		int get_no_timesteps() const;
        void set_no_timesteps(int no_timesteps);

        // ***************
        // MARK: - End condition
        // ***************
        
        /*
		int get_timepoint_zero_end_cond() const;
		void set_timepoint_zero_end_cond(int timepoint);
         */
        
        // ***************
        // MARK: - Getters (general)
        // ***************

		std::string get_name() const;
		ITptr get_ixn_param_traj() const;

        // ***************
        // MARK: - Set zero endpoint
        // ***************
        
        void set_timepoint_zero_end_cond(int timepoint);
        
        // ***************
        // MARK: - Get value at timepoint
        // ***************

		double get_val_at_timepoint(int timepoint) const;

        // ***************
        // MARK: - Solve diff eq
        // ***************

		virtual void solve_diff_eq_at_timepoint_to_minus_one(int timepoint, double dt) = 0;
        virtual void solve_diff_eq_at_timepoint_to_minus_one_l2(int timepoint, double dt, double l2_lambda, double l2_center) = 0;

        // ***************
        // MARK: - Write to file
        // ***************
        
		void write_to_file(std::string fname) const;
	};

};
