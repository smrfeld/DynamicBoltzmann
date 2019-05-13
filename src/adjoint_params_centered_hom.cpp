#include "../include/dblz_bits/adjoint_params_centered_hom.hpp"

// Other headers
#include "../include/dblz_bits/ixn_param_traj.hpp"
#include "../include/dblz_bits/ixn_param.hpp"
#include "../include/dblz_bits/general.hpp"
#include "../include/dblz_bits/diff_eq_rhs.hpp"
#include "../include/dblz_bits/moment_diff.hpp"
#include "../include/dblz_bits/center_traj.hpp"

#include <iostream>
#include <fstream>

/************************************
* Namespace for dblz
************************************/

namespace dblz {

    /****************************************
     AdjointParamsCenteredHomDerivTerm
     ****************************************/
    
    // ***************
    // MARK: - Constructor
    // ***************
    
    AdjointParamsCenteredHomDerivTerm::AdjointParamsCenteredHomDerivTerm(ITptr deriv_ixn_param_traj, std::map<int,std::map<Sptr,ITptr>> all_biases, std::map<int, std::map<Sptr, std::map<int, std::map<Sptr,ITptr>>>> all_weights, std::map<int,std::map<Sptr,CTptr>> all_center_trajs, std::map<int, std::map<int,int>> conn_mults) {
        _deriv_ixn_param_traj = deriv_ixn_param_traj;
        _all_biases = all_biases;
        _all_weights = all_weights;
        _all_center_trajs = all_center_trajs;
        _conn_mults = conn_mults;
        
        set_no_timesteps(0);
    };
    AdjointParamsCenteredHomDerivTerm::AdjointParamsCenteredHomDerivTerm(const AdjointParamsCenteredHomDerivTerm& other) {
        _copy(other);
    };
    AdjointParamsCenteredHomDerivTerm::AdjointParamsCenteredHomDerivTerm(AdjointParamsCenteredHomDerivTerm&& other) {
        _move(other);
    };
    AdjointParamsCenteredHomDerivTerm& AdjointParamsCenteredHomDerivTerm::operator=(const AdjointParamsCenteredHomDerivTerm& other) {
        if (this != &other) {
            _clean_up();
            _copy(other);
        };
        return *this;
    };
    AdjointParamsCenteredHomDerivTerm& AdjointParamsCenteredHomDerivTerm::operator=(AdjointParamsCenteredHomDerivTerm&& other) {
        if (this != &other) {
            _clean_up();
            _move(other);
        };
        return *this;
    };
    AdjointParamsCenteredHomDerivTerm::AdjointParamsCenteredHomDerivTerm()
    {
        _clean_up();
    };
    void AdjointParamsCenteredHomDerivTerm::_clean_up() {};
    void AdjointParamsCenteredHomDerivTerm::_copy(const AdjointParamsCenteredHomDerivTerm& other) {
        _deriv_ixn_param_traj = other._deriv_ixn_param_traj;
        _all_biases = other._all_biases;
        _all_weights = other._all_weights;
        _all_center_trajs = other._all_center_trajs;
        _no_timesteps = other._no_timesteps;
        _no_timepoints = other._no_timepoints;
        _vals = other._vals;
        _conn_mults = other._conn_mults;
    };
    void AdjointParamsCenteredHomDerivTerm::_move(AdjointParamsCenteredHomDerivTerm& other) {
        _deriv_ixn_param_traj = other._deriv_ixn_param_traj;
        _all_biases = other._all_biases;
        _all_weights = other._all_weights;
        _all_center_trajs = other._all_center_trajs;
        _no_timesteps = other._no_timesteps;
        _no_timepoints = other._no_timepoints;
        _vals = other._vals;
        _conn_mults = other._conn_mults;

        other._deriv_ixn_param_traj = nullptr;
        other._all_biases.clear();
        other._all_weights.clear();
        other._all_center_trajs.clear();
        other._no_timepoints = 0;
        other._no_timesteps = 0;
        other._vals.clear();
        other._conn_mults.clear();
    };

    // ***************
    // MARK: - Getters
    // ***************
    
    ITptr AdjointParamsCenteredHomDerivTerm::get_deriv_ixn_param_traj() const {
        return _deriv_ixn_param_traj;
    };
    
    // ***************
    // MARK: - Timesteps
    // ***************
    
    int AdjointParamsCenteredHomDerivTerm::get_no_timesteps() const {
        return _no_timesteps;
    };
    void AdjointParamsCenteredHomDerivTerm::set_no_timesteps(int no_timesteps) {
        _no_timesteps = no_timesteps;
        _no_timepoints = no_timesteps + 1;
        
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
    
    void AdjointParamsCenteredHomDerivTerm::calculate_val_at_timepoint(int timepoint, bool form_abscissas) {
        double val = 0.0;
        
        // Weights first term
        for (auto const &pr1: _all_weights) {
            for (auto const &pr2: pr1.second) {
                for (auto const &pr3: pr2.second) {
                    // Ensure only l-1 to l connection
                    if (pr3.first < pr1.first) {
                        continue;
                    };
                    for (auto const &pr4: pr3.second) {
                        val += pr4.second->get_adjoint()->get_val_at_timepoint(timepoint) * pr4.second->get_diff_eq_rhs()->get_deriv_wrt_nu_at_timepoint(timepoint, _deriv_ixn_param_traj,form_abscissas);
                    };
                };
            };
        };
        
        // Biases
        int layer;
        double conn_mult, center_val, adjoint_val;
        for (auto const &pr1: _all_biases) {
            layer = pr1.first;
            for (auto const &pr2: pr1.second) {
                // Adjoint val
                adjoint_val = pr2.second->get_adjoint()->get_val_at_timepoint(timepoint);
                
                // First term
                val += adjoint_val * pr2.second->get_diff_eq_rhs()->get_deriv_wrt_nu_at_timepoint(timepoint, _deriv_ixn_param_traj,form_abscissas);
                
                // Second term connecting l-1 and l
                if (layer != 0) {
                    conn_mult = _conn_mults.at(layer-1).at(layer);
                    
                    // iterate over species in l-1 layer
                    for (auto const &pr3: _all_weights.at(layer).at(pr2.first).at(layer-1)) {
                        center_val = _all_center_trajs.at(layer-1).at(pr3.first)->get_val_at_timepoint(timepoint);
                        val += adjoint_val * conn_mult * center_val * pr3.second->get_diff_eq_rhs()->get_deriv_wrt_nu_at_timepoint(timepoint, _deriv_ixn_param_traj,form_abscissas);
                    };
                };
                
                // Third term connecting l and l+1
                if (layer != _conn_mults.size()-1) {
                    conn_mult = _conn_mults.at(layer).at(layer+1);
                    
                    for (auto const &pr3: _all_weights.at(layer).at(pr2.first).at(layer+1)) {
                        center_val = _all_center_trajs.at(layer+1).at(pr3.first)->get_val_at_timepoint(timepoint);
                        val += adjoint_val * conn_mult * center_val * pr3.second->get_diff_eq_rhs()->get_deriv_wrt_nu_at_timepoint(timepoint, _deriv_ixn_param_traj,form_abscissas);
                    };
                };
            };
        };
        
        // Set
        _vals[timepoint] = val;
    };
    double AdjointParamsCenteredHomDerivTerm::get_val_at_timepoint(int timepoint) const {
        return _vals.at(timepoint);
    };

    // ***************
    // MARK: - AdjointParamsCenteredHom
    // ***************
    
    // ***************
    // MARK: - Constructor
    // ***************

    AdjointParamsCenteredHomBias::AdjointParamsCenteredHomBias(std::string name, ITptr ixn_param_traj, CtrDerivPtr deriv_term_bias) : Adjoint(name,ixn_param_traj) {
        // Check
        if (ixn_param_traj->get_type() != IxnParamType::H && ixn_param_traj->get_type() != IxnParamType::B) {
            std::cerr << ">>> AdjointParamsCenteredHomBias::AdjointParamsCenteredHomBias <<< Error: this is the wrong adjoint for weight terms." << std::endl;
            exit(EXIT_FAILURE);
        };
        
        _deriv_term_bias = deriv_term_bias;
    };
    AdjointParamsCenteredHomBias::AdjointParamsCenteredHomBias(const AdjointParamsCenteredHomBias& other) : Adjoint(other) {
		_copy(other);
	};
    AdjointParamsCenteredHomBias::AdjointParamsCenteredHomBias(AdjointParamsCenteredHomBias&& other) : Adjoint(std::move(other)) {
		_move(other);
	};
    AdjointParamsCenteredHomBias& AdjointParamsCenteredHomBias::operator=(const AdjointParamsCenteredHomBias& other) {
		if (this != &other) {
			_clean_up();
            Adjoint::operator=(other);
			_copy(other);
		};
		return *this;
    };
    AdjointParamsCenteredHomBias& AdjointParamsCenteredHomBias::operator=(AdjointParamsCenteredHomBias&& other) {
		if (this != &other) {
			_clean_up();
            Adjoint::operator=(std::move(other));
			_move(other);
		};
		return *this;
    };
	AdjointParamsCenteredHomBias::~AdjointParamsCenteredHomBias()
	{
		_clean_up();
	};
    void AdjointParamsCenteredHomBias::_clean_up() {};
	void AdjointParamsCenteredHomBias::_copy(const AdjointParamsCenteredHomBias& other) {
        _deriv_term_bias = other._deriv_term_bias;
    };
	void AdjointParamsCenteredHomBias::_move(AdjointParamsCenteredHomBias& other) {
        _deriv_term_bias = other._deriv_term_bias;
        
        other._deriv_term_bias = nullptr;
    };

    // ***************
    // MARK: - Get deriv term
    // ***************
    
    CtrDerivPtr AdjointParamsCenteredHomBias::get_deriv_term_bias() const {
        return _deriv_term_bias;
    };

    // ***************
    // MARK: - Diff eq
    // ***************
    
	void AdjointParamsCenteredHomBias::solve_diff_eq_at_timepoint_to_minus_one(int timepoint, double dt, bool form_abscissas) {
        
		// Deriv term
		double deriv_term = _deriv_term_bias->get_val_at_timepoint(timepoint);
        
		// Difference in moments
        double moment_delta = -1.0 * _ixn_param_traj->get_ixn_param_at_timepoint(timepoint)->get_moment_diff()->get_moment_diff_awake_minus_asleep_plus_offset();
        
        // Step
        _vals[timepoint-1] = _vals[timepoint] - dt * (moment_delta - deriv_term);
	};
    
    void AdjointParamsCenteredHomBias::solve_diff_eq_at_timepoint_to_minus_one_l2(int timepoint, double dt, double l2_lambda, double l2_center, bool form_abscissas) {
        
        std::cerr << ">>> AdjointParamsCenteredHomBias::solve_diff_eq_at_timepoint_to_minus_one_l2 <<< Error: L2 mode not yet supported" << std::endl;
        exit(EXIT_FAILURE);
    };
    
    // ***************
    // MARK: - AdjointParamsCenteredHomWeight
    // ***************
    
    // ***************
    // MARK: - Constructor
    // ***************
    
    AdjointParamsCenteredHomWeight::AdjointParamsCenteredHomWeight(std::string name, ITptr ixn_param_traj, CtrDerivPtr deriv_term_weight, CtrDerivPtr deriv_term_bias_lower, CtrDerivPtr deriv_term_bias_upper, int conn_mult, CTptr center_lower, CTptr center_upper, std::shared_ptr<AdjointParamsCenteredHomBias> adjoint_bias_lower, std::shared_ptr<AdjointParamsCenteredHomBias> adjoint_bias_upper, std::vector<CTptr> all_centers_lower, std::vector<CTptr> all_centers_upper) : Adjoint(name,ixn_param_traj) {
        if (ixn_param_traj->get_type() != IxnParamType::W && ixn_param_traj->get_type() != IxnParamType::X) {
            std::cerr << ">>> AdjointParamsCenteredHomWeight::AdjointParamsCenteredHomWeight <<< Error: this is the wrong constructor for bias terms." << std::endl;
            exit(EXIT_FAILURE);
        };
        
        _deriv_term_weight = deriv_term_weight;
        _deriv_term_bias_lower = deriv_term_bias_lower;
        _deriv_term_bias_upper = deriv_term_bias_upper;
        _conn_mult = conn_mult;
        _center_lower = center_lower;
        _center_upper = center_upper;
        _all_centers_lower = all_centers_lower;
        _all_centers_upper = all_centers_upper;
        _adjoint_bias_lower = adjoint_bias_lower;
        _adjoint_bias_upper = adjoint_bias_upper;
    };
    AdjointParamsCenteredHomWeight::AdjointParamsCenteredHomWeight(const AdjointParamsCenteredHomWeight& other) : Adjoint(other) {
        _copy(other);
    };
    AdjointParamsCenteredHomWeight::AdjointParamsCenteredHomWeight(AdjointParamsCenteredHomWeight&& other) : Adjoint(std::move(other)) {
        _move(other);
    };
    AdjointParamsCenteredHomWeight& AdjointParamsCenteredHomWeight::operator=(const AdjointParamsCenteredHomWeight& other) {
        if (this != &other) {
            _clean_up();
            Adjoint::operator=(other);
            _copy(other);
        };
        return *this;
    };
    AdjointParamsCenteredHomWeight& AdjointParamsCenteredHomWeight::operator=(AdjointParamsCenteredHomWeight&& other) {
        if (this != &other) {
            _clean_up();
            Adjoint::operator=(std::move(other));
            _move(other);
        };
        return *this;
    };
    AdjointParamsCenteredHomWeight::~AdjointParamsCenteredHomWeight()
    {
        _clean_up();
    };
    void AdjointParamsCenteredHomWeight::_clean_up() {};
    void AdjointParamsCenteredHomWeight::_copy(const AdjointParamsCenteredHomWeight& other) {
        _deriv_term_weight = other._deriv_term_weight;
        _deriv_term_bias_lower = other._deriv_term_bias_lower;
        _deriv_term_bias_upper = other._deriv_term_bias_upper;
        _conn_mult = other._conn_mult;
        _center_lower = other._center_lower;
        _center_upper = other._center_upper;
        _all_centers_lower = other._all_centers_lower;
        _all_centers_upper = other._all_centers_upper;
        _adjoint_bias_lower = other._adjoint_bias_lower;
        _adjoint_bias_upper = other._adjoint_bias_upper;
    };
    void AdjointParamsCenteredHomWeight::_move(AdjointParamsCenteredHomWeight& other) {
        _deriv_term_weight = other._deriv_term_weight;
        _deriv_term_bias_lower = other._deriv_term_bias_lower;
        _deriv_term_bias_upper = other._deriv_term_bias_upper;
        _conn_mult = other._conn_mult;
        _center_lower = other._center_lower;
        _center_upper = other._center_upper;
        _all_centers_lower = other._all_centers_lower;
        _all_centers_upper = other._all_centers_upper;
        _adjoint_bias_lower = other._adjoint_bias_lower;
        _adjoint_bias_upper = other._adjoint_bias_upper;

        other._deriv_term_weight = nullptr;
        other._deriv_term_bias_lower = nullptr;
        other._deriv_term_bias_upper = nullptr;
        other._conn_mult = 0;
        other._center_lower = nullptr;
        other._center_upper = nullptr;
        other._all_centers_lower.clear();
        other._all_centers_upper.clear();
        other._adjoint_bias_lower = nullptr;
        other._adjoint_bias_upper = nullptr;
    };
    
    // ***************
    // MARK: - Get deriv term
    // ***************
    
    CtrDerivPtr AdjointParamsCenteredHomWeight::get_deriv_term_weight() const {
        return _deriv_term_weight;
    };
    CtrDerivPtr AdjointParamsCenteredHomWeight::get_deriv_term_bias_lower() const {
        return _deriv_term_bias_lower;
    };
    CtrDerivPtr AdjointParamsCenteredHomWeight::get_deriv_term_bias_upper() const {
        return _deriv_term_bias_upper;
    };
    
    // ***************
    // MARK: - Solve diff eq
    // ***************
    
    void AdjointParamsCenteredHomWeight::solve_diff_eq_at_timepoint_to_minus_one(int timepoint, double dt, bool form_abscissas) {
        
        // Difference in moments
        double moment_delta = -1.0 * _ixn_param_traj->get_ixn_param_at_timepoint(timepoint)->get_moment_diff()->get_moment_diff_awake_minus_asleep_plus_offset();
        
        // Derivative in centers terms
        double deriv_centers_terms = 0.0;
        double f;
        
        // Lower - upper
        f = _adjoint_bias_lower->get_val_at_timepoint(timepoint) * _conn_mult;
        for (auto center_upper: _all_centers_upper) {
            deriv_centers_terms += f * center_upper->get_deriv_at_timepoint(timepoint, dt);
        };
        // Upper -lower
        f = _adjoint_bias_upper->get_val_at_timepoint(timepoint) * _conn_mult;
        for (auto center_lower: _all_centers_lower) {
            deriv_centers_terms += f * center_lower->get_deriv_at_timepoint(timepoint, dt);
        };
        
        // Deriv term weight
        double deriv_term = _deriv_term_weight->get_val_at_timepoint(timepoint);
        
        // Other deriv terms
        double deriv_mixing_terms = _conn_mult * _center_upper->get_val_at_timepoint(timepoint) * _deriv_term_bias_lower->get_val_at_timepoint(timepoint) + _conn_mult * _center_lower->get_val_at_timepoint(timepoint) * _deriv_term_bias_upper->get_val_at_timepoint(timepoint);
        
        // Step
        _vals[timepoint-1] = _vals[timepoint] - dt * (moment_delta - deriv_centers_terms - deriv_term + deriv_mixing_terms);
    };
    void AdjointParamsCenteredHomWeight::solve_diff_eq_at_timepoint_to_minus_one_l2(int timepoint, double dt, double l2_lambda, double l2_center, bool form_abscissas) {
        
        std::cerr << ">>> AdjointParamsCenteredHomWeight::solve_diff_eq_at_timepoint_to_minus_one_l2 <<< Error: L2 mode not yet supported" << std::endl;
        exit(EXIT_FAILURE);
    };
};
