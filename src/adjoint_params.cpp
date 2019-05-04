#include "../include/dblz_bits/adjoint_params.hpp"

// Other headers
#include "../include/dblz_bits/ixn_param_traj.hpp"
#include "../include/dblz_bits/ixn_param.hpp"
#include "../include/dblz_bits/general.hpp"
#include "../include/dblz_bits/diff_eq_rhs.hpp"
#include "../include/dblz_bits/moment_diff.hpp"

#include <iostream>
#include <fstream>

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	AdjointParams
	****************************************/

    // ***************
    // MARK: - Constructor
    // ***************

    AdjointParams::AdjointParams(std::string name, ITptr ixn_param_traj) : Adjoint(name,ixn_param_traj) {};
    AdjointParams::AdjointParams(const AdjointParams& other) : Adjoint(other) {
		_copy(other);
	};
    AdjointParams::AdjointParams(AdjointParams&& other) : Adjoint(std::move(other)) {
		_move(other);
	};
    AdjointParams& AdjointParams::operator=(const AdjointParams& other) {
		if (this != &other) {
			_clean_up();
            Adjoint::operator=(other);
			_copy(other);
		};
		return *this;
    };
    AdjointParams& AdjointParams::operator=(AdjointParams&& other) {
		if (this != &other) {
			_clean_up();
            Adjoint::operator=(std::move(other));
			_move(other);
		};
		return *this;
    };
	AdjointParams::~AdjointParams()
	{
		_clean_up();
	};
    void AdjointParams::_clean_up() {};
	void AdjointParams::_copy(const AdjointParams& other) {};
	void AdjointParams::_move(AdjointParams& other) {};

    // ***************
    // MARK: - Diff eq
    // ***************
    
	void AdjointParams::solve_diff_eq_at_timepoint_to_minus_one(int timepoint, double dt) {
		if (timepoint >= _no_timepoints) {
			std::cerr << ">>> Error: AdjointParams::solve_diff_eq_at_timepoint_to_minus_one <<< " << timepoint << " is out of bounds: " << _no_timepoints << std::endl;
			exit(EXIT_FAILURE);
		};
		if (timepoint > _timepoint_zero_end_cond) {
			std::cerr << ">>> Error: AdjointParams::solve_diff_eq_at_timepoint_to_minus_one <<< " << timepoint << " is beyond the zero endpoint: " << _timepoint_zero_end_cond << std::endl;
			exit(EXIT_FAILURE);		
		};

		// Deriv val
		double deriv_term=0.0;
		double deriv, adjoint_val;

		// Go through dependencies (numerator F and AdjointParams)
        for (auto const dep_pair: _ixn_param_traj->get_diff_eq_dependencies()) {

			// Derivative
			deriv = dep_pair.first->get_deriv_wrt_nu_at_timepoint(timepoint,dep_pair.second);

			// AdjointParams
			auto adjoint = dep_pair.first->get_parent_ixn_param_traj()->get_adjoint();
			if (!adjoint) {
				std::cerr << ">>> Error: AdjointParams::solve_diff_eq_at_timepoint_to_minus_one <<< No AdjointParams for ixn param: " << dep_pair.first->get_parent_ixn_param_traj()->get_name() << std::endl;
				exit(EXIT_FAILURE);	
			};
			adjoint_val = adjoint->get_val_at_timepoint(timepoint);

			// Add
			deriv_term += deriv * adjoint_val;

		};

		// Difference in moments
        auto moment = _ixn_param_traj->get_ixn_param_at_timepoint(timepoint)->get_moment();
        double moment_delta = -1.0 * moment->get_moment_diff_awake_minus_asleep_plus_offset();
        
        // Step
        _vals[timepoint-1] = _vals[timepoint] - dt * (moment_delta - deriv);
	};
    
    void AdjointParams::solve_diff_eq_at_timepoint_to_minus_one_l2(int timepoint, double dt, double l2_lambda, double l2_center) {

        // Solve
        solve_diff_eq_at_timepoint_to_minus_one(timepoint, dt);
        
        // L2 reg
        double ixn_param_val = _ixn_param_traj->get_ixn_param_at_timepoint(timepoint)->get_val();
        double l2_term = 2.0 * l2_lambda * (ixn_param_val - l2_center);
        
        // Step
        _vals[timepoint-1] -= dt * l2_term;
    };
};
