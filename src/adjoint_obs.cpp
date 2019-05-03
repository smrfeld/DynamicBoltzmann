#include "../include/dblz_bits/adjoint_obs.hpp"

// Other headers
#include "../include/dblz_bits/ixn_param_traj.hpp"
#include "../include/dblz_bits/ixn_param.hpp"
#include "../include/dblz_bits/general.hpp"
#include "../include/dblz_bits/diff_eq_rhs.hpp"
#include "../include/dblz_bits/moment.hpp"

#include <iostream>
#include <fstream>

/************************************
* Namespace for dblz
************************************/

namespace dblz {
    
    // ***************
    // MARK: - AdjointObs
    // ***************

    // ***************
    // MARK: - Constructor
    // ***************

    AdjointObs::AdjointObs(std::string name, ITptr ixn_param_traj) : Adjoint(name, ixn_param_traj) {};
    AdjointObs::AdjointObs(const AdjointObs& other) : Adjoint(other) {
		_copy(other);
	};
    AdjointObs::AdjointObs(AdjointObs&& other) : Adjoint(std::move(other)) {
		_move(other);
	};
    AdjointObs& AdjointObs::operator=(const AdjointObs& other) {
		if (this != &other) {
			_clean_up();
            Adjoint::operator=(other);
			_copy(other);
		};
		return *this;
    };
    AdjointObs& AdjointObs::operator=(AdjointObs&& other) {
		if (this != &other) {
			_clean_up();
            Adjoint::operator=(std::move(other));
			_move(other);
		};
		return *this;
    };
	AdjointObs::~AdjointObs()
	{
		_clean_up();
	};
	void AdjointObs::_clean_up() {};
	void AdjointObs::_copy(const AdjointObs& other) {};
	void AdjointObs::_move(AdjointObs& other) {};

    // ***************
    // MARK: - Diff eq
    // ***************
    
	void AdjointObs::solve_diff_eq_at_timepoint_to_minus_one(int timepoint, double dt) {
		if (timepoint >= _no_timepoints) {
			std::cerr << ">>> Error: AdjointObs::solve_diff_eq_at_timepoint_to_minus_one <<< " << timepoint << " is out of bounds: " << _no_timepoints << std::endl;
			exit(EXIT_FAILURE);
		};
		if (timepoint > _timepoint_zero_end_cond) {
			std::cerr << ">>> Error: AdjointObs::solve_diff_eq_at_timepoint_to_minus_one <<< " << timepoint << " is beyond the zero endpoint: " << _timepoint_zero_end_cond << std::endl;
			exit(EXIT_FAILURE);		
		};

		// Deriv val
		double deriv_term=0.0;
		double deriv, adjoint_val;

		// Go through dependencies (numerator F and AdjointObs)
        for (auto const dep_pair: _ixn_param_traj->get_diff_eq_dependencies()) {

			// Derivative
			deriv = dep_pair.first->get_deriv_wrt_nu_at_timepoint(timepoint,dep_pair.second);

			// AdjointObs
			auto adjoint = dep_pair.first->get_parent_ixn_param_traj()->get_adjoint();
			if (!adjoint) {
				std::cerr << ">>> Error: AdjointObs::solve_diff_eq_at_timepoint_to_minus_one <<< No AdjointObs for ixn param: " << dep_pair.first->get_parent_ixn_param_traj()->get_name() << std::endl;
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
    
    void AdjointObs::solve_diff_eq_at_timepoint_to_minus_one_l2(int timepoint, double dt, double l2_lambda, double l2_center) {

        // Solve
        solve_diff_eq_at_timepoint_to_minus_one(timepoint, dt);
        
        // L2 reg
        double ixn_param_val = _ixn_param_traj->get_ixn_param_at_timepoint(timepoint)->get_val();
        double l2_term = 2.0 * l2_lambda * (ixn_param_val - l2_center);
        
        // Step
        _vals[timepoint-1] -= dt * l2_term;
    };
    
    // ***************
    // MARK: - AdjointObsCommonTerm
    // ***************
    
    // ***************
    // MARK: - Constructor
    // ***************
    
    AdjointObsCommonTerm::AdjointObsCommonTerm(Sptr species, int idx_diff_eq, std::vector<ITptr> all_ixn_param_trajs) {
        _species = species;
        _idx_diff_eq = idx_diff_eq;
        _all_ixn_param_trajs = all_ixn_param_trajs;
        
        set_no_timesteps(0);
    };
    AdjointObsCommonTerm::AdjointObsCommonTerm(const AdjointObsCommonTerm& other) {
        _copy(other);
    };
    AdjointObsCommonTerm::AdjointObsCommonTerm(AdjointObsCommonTerm&& other) {
        _move(other);
    };
    AdjointObsCommonTerm& AdjointObsCommonTerm::operator=(const AdjointObsCommonTerm& other) {
        if (this != &other) {
            _clean_up();
            _copy(other);
        };
        return *this;
    };
    AdjointObsCommonTerm& AdjointObsCommonTerm::operator=(AdjointObsCommonTerm&& other) {
        if (this != &other) {
            _clean_up();
            _move(other);
        };
        return *this;
    };
    AdjointObsCommonTerm::~AdjointObsCommonTerm()
    {
        _clean_up();
    };
    void AdjointObsCommonTerm::_clean_up() {};
    void AdjointObsCommonTerm::_copy(const AdjointObsCommonTerm& other) {
        _species = other._species;
        _idx_diff_eq = other._idx_diff_eq;
        _all_ixn_param_trajs = other._all_ixn_param_trajs;
        _no_timesteps = other._no_timesteps;
        _no_timepoints = other._no_timepoints;
        _vals = other._vals;
    };
    void AdjointObsCommonTerm::_move(AdjointObsCommonTerm& other) {
        _species = other._species;
        _idx_diff_eq = other._idx_diff_eq;
        _all_ixn_param_trajs = other._all_ixn_param_trajs;
        _no_timesteps = other._no_timesteps;
        _no_timepoints = other._no_timepoints;
        _vals = other._vals;

        other._species = nullptr;
        other._idx_diff_eq = 0;
        other._all_ixn_param_trajs.clear();
        other._no_timepoints = 0;
        other._no_timesteps = 0;
        other._vals.clear();
    };
    
    // ***************
    // MARK: - Timesteps
    // ***************
    
    int AdjointObsCommonTerm::get_no_timesteps() const {
        return _no_timesteps;
    };
    void AdjointObsCommonTerm::set_no_timesteps(int no_timesteps) {
        _no_timesteps = no_timesteps;
        _no_timepoints = no_timesteps+1;
        
        while (_vals.size() < _no_timepoints) {
            _vals.push_back(0.0);
        };
        while (_vals.size() > _no_timepoints) {
            _vals.pop_back();
        };
    };
    
    // ***************
    // MARK: - Vals
    // ***************
    
    void AdjointObsCommonTerm::calculate_val_at_timepoint(int timepoint) {
        _vals[timepoint] = 0.0;
        
        for (auto ixn: _all_ixn_param_trajs) {
            auto adjoint = ixn->get_adjoint();
            auto diff_eq_rhs = ixn->get_diff_eq_rhs();
            
            _vals[timepoint] += ixn->get_adjoint()->get_val_at_timepoint(timepoint) * ixn->get_diff_eq_rhs()->get_deriv_wrt_nu_at_timepoint(timepoint, _idx_diff_eq);
        };
    };
    double AdjointObsCommonTerm::get_val_at_timepoint(int timepoint) const {
        return _vals.at(timepoint);
    };
};
