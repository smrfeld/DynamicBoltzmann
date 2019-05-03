#include <string>
#include <map>
#include <vector>

#include "fwds/fwds_ixn_param_traj.hpp"

#ifndef ADJOINT_H
#define ADJOINT_H
#include "adjoint.hpp"
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {	

    // ***************
    // MARK: - Adjoint class when params used for diff eq RHS
    // ***************
    
    class AdjointParams : public Adjoint {

	private:

		// Internal copy func/clean up
		void _clean_up();
		void _copy(const AdjointParams& other);
		void _move(AdjointParams &other);

	public:

        // ***************
        // MARK: - Constructor
        // ***************

		AdjointParams(std::string name, ITptr ixn_param_traj);
		AdjointParams(const AdjointParams& other);
		AdjointParams& operator=(const AdjointParams& other);
		AdjointParams(AdjointParams&& other);
		AdjointParams& operator=(AdjointParams&& other);
		~AdjointParams();

        // ***************
        // MARK: - Solve diff eq
        // ***************

		void solve_diff_eq_at_timepoint_to_minus_one(int timepoint, double dt);
        void solve_diff_eq_at_timepoint_to_minus_one_l2(int timepoint, double dt, double l2_lambda, double l2_center);
	};

};
