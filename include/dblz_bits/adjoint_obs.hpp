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
    // MARK: - Adjoint class when observables used for diff eq RHS
    // ***************
    
    class AdjointObs : public Adjoint {

	private:

		// Internal copy func/clean up
		void _clean_up();
		void _copy(const AdjointObs& other);
		void _move(AdjointObs &other);

	public:

        // ***************
        // MARK: - Constructor
        // ***************

		AdjointObs(std::string name, ITptr ixn_param_traj);
		AdjointObs(const AdjointObs& other);
		AdjointObs& operator=(const AdjointObs& other);
		AdjointObs(AdjointObs&& other);
		AdjointObs& operator=(AdjointObs&& other);
		~AdjointObs();

        // ***************
        // MARK: - Solve diff eq
        // ***************

        void solve_diff_eq_at_timepoint_to_minus_one(int timepoint, double dt);
        void solve_diff_eq_at_timepoint_to_minus_one_l2(int timepoint, double dt, double l2_lambda, double l2_center);
	};

};
