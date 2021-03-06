#include "../include/dblz_bits/ixn_param_traj.hpp"

// Other headers
#include "../include/dblz_bits/general.hpp"
#include "../include/dblz_bits/species.hpp"
#include "../include/dblz_bits/diff_eq_rhs.hpp"
#include "../include/dblz_bits/adjoint_obs.hpp"
#include "../include/dblz_bits/adjoint_params.hpp"
#include "../include/dblz_bits/adjoint_params_centered_hom.hpp"
#include "../include/dblz_bits/ixn_param.hpp"
#include "../include/dblz_bits/moment_diff.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	IxnFunc
	****************************************/

    // ***************
    // MARK: - Constructor
    // ***************

	IxnParamTraj::IxnParamTraj(std::string name, IxnParamType type, double init_cond) {
		_init_cond = init_cond;

        // First element
        _ixn_params.push_back(std::make_shared<IxnParam>(name,type,init_cond,0.0));
        
        // Adjoint
        _adjoint = nullptr;
        _adjoint_obs = nullptr;
        _adjoint_params = nullptr;
        _adjoint_params_centered_hom_bias = nullptr;
        _adjoint_params_centered_hom_weight = nullptr;

        // Diff eq
        _diff_eq = nullptr;
        
        // No timesteps
        set_no_timesteps(0);
        
        _is_val_fixed = false;
	};
	IxnParamTraj::IxnParamTraj(const IxnParamTraj& other) {
		_copy(other);
	};
	IxnParamTraj::IxnParamTraj(IxnParamTraj&& other) {
		_move(other);
	};
    IxnParamTraj& IxnParamTraj::operator=(const IxnParamTraj& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    IxnParamTraj& IxnParamTraj::operator=(IxnParamTraj&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };
	IxnParamTraj::~IxnParamTraj()
	{
		_clean_up();
	};
	void IxnParamTraj::_clean_up() {
        // Nothing...
    };
	void IxnParamTraj::_copy(const IxnParamTraj& other) {
        _adjoint = other._adjoint;
        _adjoint_obs = other._adjoint_obs;
        _adjoint_params = other._adjoint_params;
        _adjoint_params_centered_hom_bias = other._adjoint_params_centered_hom_bias;
        _adjoint_params_centered_hom_weight = other._adjoint_params_centered_hom_weight;
        
        _diff_eq = other._diff_eq;
        
        _diff_eq_dependencies = other._diff_eq_dependencies;
        
        _ixn_params = other._ixn_params;
        
        _substep_val_current = other._substep_val_current;
        _substep_val_new = other._substep_val_new;
        
        _init_cond = other._init_cond;
        
        _no_timesteps = other._no_timesteps;
        _no_timepoints = other._no_timepoints;
        
        _is_val_fixed = other._is_val_fixed;
    };
	void IxnParamTraj::_move(IxnParamTraj& other) {
        _adjoint = other._adjoint;
        _adjoint_obs = other._adjoint_obs;
        _adjoint_params = other._adjoint_params;
        _adjoint_params_centered_hom_bias = other._adjoint_params_centered_hom_bias;
        _adjoint_params_centered_hom_weight = other._adjoint_params_centered_hom_weight;

        _diff_eq = other._diff_eq;
        
        _diff_eq_dependencies = other._diff_eq_dependencies;
        
        _ixn_params = other._ixn_params;
        
        _substep_val_current = other._substep_val_current;
        _substep_val_new = other._substep_val_new;

        _init_cond = other._init_cond;
        
        _no_timesteps = other._no_timesteps;
        _no_timepoints = other._no_timepoints;
        
        _is_val_fixed = other._is_val_fixed;

		// Reset the other
        other._adjoint = nullptr;
        other._adjoint_obs = nullptr;
        other._adjoint_params = nullptr;
        other._adjoint_params_centered_hom_bias = nullptr;
        other._adjoint_params_centered_hom_weight = nullptr;
        other._diff_eq = nullptr;
        other._diff_eq_dependencies.clear();
        other._ixn_params.clear();
        other._substep_val_current = 0.0;
        other._substep_val_new = 0.0;
        other._init_cond = 0.0;
        other._no_timesteps = 0;
        other._no_timepoints = 1;
        other._is_val_fixed = false;
    };

    // ***************
    // MARK: - Print moment string
    // ***************
    
    void IxnParamTraj::print_val_traj(int timepoint_start, int no_timesteps) const {
        for (auto timepoint=timepoint_start; timepoint<=timepoint_start+no_timesteps; timepoint++) {
            std::cout << _ixn_params.at(timepoint)->get_val();
            if (timepoint != timepoint_start + no_timesteps) {
                std::cout << " ";
            };
        };
        std::cout << std::endl;
    };
    
    void IxnParamTraj::print_moment_traj(int timepoint_start, int no_timesteps) const {
        for (auto timepoint=timepoint_start; timepoint<=timepoint_start+no_timesteps; timepoint++) {
            std::cout << _ixn_params.at(timepoint)->get_moment_diff()->get_moment_comparison_str();
            if (timepoint != timepoint_start+no_timesteps) {
                std::cout << " ";
            };
        };
        std::cout << std::endl;
    };
    
    void IxnParamTraj::print_moment_diff_traj(int timepoint_start, int no_timesteps) const {
        for (auto timepoint=timepoint_start; timepoint<=timepoint_start+no_timesteps; timepoint++) {
            std::cout << _ixn_params.at(timepoint)->get_moment_diff()->get_moment(MCType::AWAKE);
            // std::cout << _diff_eq->get_lr() * ( _ixn_params.at(timepoint)->get_moment_diff()->get_moment(MCType::AWAKE) - _ixn_params.at(timepoint)->get_moment_diff()->get_moment(MCType::ASLEEP) );
            if (timepoint != timepoint_start+no_timesteps) {
                std::cout << " ";
            };
        };
        std::cout << std::endl;

        for (auto timepoint=timepoint_start; timepoint<=timepoint_start+no_timesteps; timepoint++) {
            std::cout << _ixn_params.at(timepoint)->get_moment_diff()->get_moment(MCType::ASLEEP);
            if (timepoint != timepoint_start+no_timesteps) {
                std::cout << " ";
            };
        };
        std::cout << std::endl;
    };
    
    // ***************
    // MARK: - Substep vals
    // ***************
    
    double IxnParamTraj::get_substep_val() const {
        return _substep_val_current;
    };
    void IxnParamTraj::set_substep_val_new(double val) {
        _substep_val_new = val;
    };
    void IxnParamTraj::set_substep_val_current(double val) {
        _substep_val_current = val;
    };
    void IxnParamTraj::move_to_new_substep_val() {
        _substep_val_current = _substep_val_new;
    };
    
    // ***************
    // MARK: - Diff eq rhs that this ixn param appears in
    // ***************
    
	void IxnParamTraj::add_diff_eq_dependency(std::shared_ptr<DiffEqRHS> diff_eq, int idx) {
		_diff_eq_dependencies.push_back(std::make_pair(diff_eq,idx));
	};
	const std::vector<std::pair<std::shared_ptr<DiffEqRHS>,int>>& IxnParamTraj::get_diff_eq_dependencies() const {
		return _diff_eq_dependencies;
	};
    
    // ***************
    // MARK: - Fix
    // ***************
    
    void IxnParamTraj::set_fix_value(bool fixed) {
        _is_val_fixed = fixed;
    };
    bool IxnParamTraj::get_is_val_fixed() const {
        return _is_val_fixed;
    };

    // ***************
    // MARK: - Timesteps
    // ***************
    
	int IxnParamTraj::get_no_timesteps() const {
		return _no_timesteps;
	};
	void IxnParamTraj::set_no_timesteps(int no_timesteps) {
		_no_timesteps = no_timesteps;
		_no_timepoints = _no_timesteps + 1;

        if (_no_timepoints <= 0) {
            std::cerr << ">>> IxnParamTraj::set_no_timesteps <<< Error: no_timepoints must be > 0. Currently, it is: " << _no_timepoints << std::endl;
            exit(EXIT_FAILURE);
        };
        
		// Adjust ixn params
        while (_ixn_params.size() < _no_timepoints) {
            // Add
            _ixn_params.push_back(std::make_shared<IxnParam>(*_ixn_params.back()));
        };
        while (_ixn_params.size() > _no_timepoints) {
            // Remove
            _ixn_params.pop_back();
        };
        
        // Set for adjoints of ixn params
        if (_adjoint) {
            _adjoint->set_no_timesteps(no_timesteps);
        };
        if (_adjoint_params_centered_hom_bias) {
            if (_adjoint_params_centered_hom_bias->get_deriv_term_bias()) { _adjoint_params_centered_hom_bias->get_deriv_term_bias()->set_no_timesteps(no_timesteps);
            };
        };
        if (_adjoint_params_centered_hom_weight) {
            if (_adjoint_params_centered_hom_weight->get_deriv_term_weight()) {
                _adjoint_params_centered_hom_weight->get_deriv_term_weight()->set_no_timesteps(no_timesteps);
            };
           
            if (_adjoint_params_centered_hom_weight->get_deriv_term_bias_upper()) {
               _adjoint_params_centered_hom_weight->get_deriv_term_bias_upper()->set_no_timesteps(no_timesteps);
            };
            
            if (_adjoint_params_centered_hom_weight->get_deriv_term_bias_lower()) {
                _adjoint_params_centered_hom_weight->get_deriv_term_bias_lower()->set_no_timesteps(no_timesteps);
            };
        };
        
        // Set for domain of diff eq rhs if needed
        if (_diff_eq) {
            for (auto domain_obs: _diff_eq->get_domain()->get_domain_obs()) {
                domain_obs->set_no_timesteps(no_timesteps);
            };
        };
        
        // Set for common/covariance terms in adjoint obs
        if (_adjoint_obs) {
            for (auto term_pr: _adjoint_obs->get_terms()) {
                term_pr.first->set_no_timesteps(no_timesteps);
                term_pr.second->set_no_timesteps(no_timesteps);
            };
        };
    };
    
    // ***************
    // MARK: - Init cond
    // ***************
    
	double IxnParamTraj::get_init_cond() const {
		return _init_cond;
	};
	void IxnParamTraj::set_init_cond(double init_cond) {
		_init_cond = init_cond;
        if (_no_timepoints >= 1) {
            _ixn_params.front()->set_val(_init_cond);
        };
    };

    // ***************
    // MARK: - Name, type
    // ***************

	std::string IxnParamTraj::get_name() const {
        return _ixn_params.front()->get_name();
	};

	IxnParamType IxnParamTraj::get_type() const {
        return _ixn_params.front()->get_type();
	};

    // ***************
    // MARK: - Diff eq
    // ***************
    
	std::shared_ptr<DiffEqRHS> IxnParamTraj::get_diff_eq_rhs() const {
		return _diff_eq;
	};
	void IxnParamTraj::set_diff_eq_rhs(std::shared_ptr<DiffEqRHS> diff_eq) {
		_diff_eq = diff_eq;
	};

	void IxnParamTraj::solve_diff_eq_at_timepoint_to_plus_one(int timepoint, double dt, bool form_abscissas) {
        /*
		if (timepoint >= _no_timepoints) {
			std::cerr << ">>> Error: IxnParamTraj::solve_diff_eq_at_timepoint_to_plus_one <<< " << timepoint << " is out of bounds: " << _no_timepoints << std::endl;
			exit(EXIT_FAILURE);
		};
        
		if (!_diff_eq) {
			std::cerr << ">>> Error: IxnParamTraj::solve_diff_eq_at_timepoint_to_plus_one <<< No diff eq set" << std::endl;
			exit(EXIT_FAILURE);
		};
         */

        _ixn_params.at(timepoint+1)->set_val(_ixn_params.at(timepoint)->get_val() + dt * _diff_eq->get_val_at_timepoint(timepoint,form_abscissas));
	};

    // ***************
    // MARK: - Ixn param
    // ***************
    
    std::shared_ptr<IxnParam> IxnParamTraj::get_ixn_param_at_timepoint(int timepoint) const {
        return _ixn_params.at(timepoint);
    };
    
    // ***************
    // MARK: - Adjoint
    // ***************
    
	void IxnParamTraj::set_adjoint(std::shared_ptr<AdjointObs> adjoint_obs) {
		_adjoint_obs = adjoint_obs;
        _adjoint = _adjoint_obs;
        _adjoint_params = nullptr;
        _adjoint_params_centered_hom_bias = nullptr;
        _adjoint_params_centered_hom_weight = nullptr;
	};
    void IxnParamTraj::set_adjoint(std::shared_ptr<AdjointParams> adjoint_params) {
        _adjoint_obs = nullptr;
        _adjoint = adjoint_params;
        _adjoint_params = adjoint_params;
        _adjoint_params_centered_hom_bias = nullptr;
        _adjoint_params_centered_hom_weight = nullptr;
    };
    void IxnParamTraj::set_adjoint(std::shared_ptr<AdjointParamsCenteredHomBias> adjoint_params_centered_hom_bias) {
        _adjoint_obs = nullptr;
        _adjoint = adjoint_params_centered_hom_bias;
        _adjoint_params = nullptr;
        _adjoint_params_centered_hom_bias = adjoint_params_centered_hom_bias;
        _adjoint_params_centered_hom_weight = nullptr;
    };
    void IxnParamTraj::set_adjoint(std::shared_ptr<AdjointParamsCenteredHomWeight> adjoint_params_centered_hom_weight) {
        _adjoint_obs = nullptr;
        _adjoint = adjoint_params_centered_hom_weight;
        _adjoint_params = nullptr;
        _adjoint_params_centered_hom_bias = nullptr;
        _adjoint_params_centered_hom_weight = adjoint_params_centered_hom_weight;
    };
	std::shared_ptr<Adjoint> IxnParamTraj::get_adjoint() const {
		return _adjoint;
	};
    std::shared_ptr<AdjointObs> IxnParamTraj::get_adjoint_obs() const {
        return _adjoint_obs;
    };
    std::shared_ptr<AdjointParams> IxnParamTraj::get_adjoint_params() const {
        return _adjoint_params;
    };
    std::shared_ptr<AdjointParamsCenteredHomBias> IxnParamTraj::get_adjoint_params_centered_hom_bias() const {
        return _adjoint_params_centered_hom_bias;
    };
    std::shared_ptr<AdjointParamsCenteredHomWeight> IxnParamTraj::get_adjoint_params_centered_hom_weight() const {
        return _adjoint_params_centered_hom_weight;
    };

    // ***************
    // MARK: - Write to file
    // ***************

	void IxnParamTraj::write_val_traj_to_file(int timepoint_start, int no_timesteps, std::string fname, bool with_timepoints) const {
		std::ofstream f;

		// Open
		f.open(fname);

		// Make sure we found it
		if (!f.is_open()) {
			std::cerr << ">>> Error: IxnParamTraj::write_val_traj_to_file <<< could not write to file: " << fname << std::endl;
			exit(EXIT_FAILURE);
		};

		// Go through all time
        if (with_timepoints) {
            for (auto timepoint=timepoint_start; timepoint<=timepoint_start+no_timesteps; timepoint++) {
                f << timepoint << " " << _ixn_params.at(timepoint)->get_val() << "\n";
            };
        } else {
            for (auto timepoint=timepoint_start; timepoint<=timepoint_start+no_timesteps; timepoint++) {
                f << _ixn_params.at(timepoint)->get_val() << "\n";
            };
        };
        
		// Close
		f.close();
	};
    
    void IxnParamTraj::write_moment_traj_to_file(int timepoint_start, int no_timesteps, std::string fname, bool with_timepoints) const {
        std::ofstream f;
        
        // Open
        f.open(fname);
        
        // Make sure we found it
        if (!f.is_open()) {
            std::cerr << ">>> Error: IxnParamTraj::write_moment_traj_to_file <<< could not write to file: " << fname << std::endl;
            exit(EXIT_FAILURE);
        };
        
        // Go through all time
        std::shared_ptr<MomentDiff> moment;
        if (with_timepoints) {
            for (auto timepoint=timepoint_start; timepoint<=timepoint_start+no_timesteps; timepoint++) {
                moment = _ixn_params.at(timepoint)->get_moment_diff();
                f << timepoint << " " << moment->get_moment(MCType::AWAKE) << " " << moment->get_moment(MCType::ASLEEP) << "\n";
            };
        } else {
            for (auto timepoint=timepoint_start; timepoint<=timepoint_start+no_timesteps; timepoint++) {
                moment = _ixn_params.at(timepoint)->get_moment_diff();
                f << moment->get_moment(MCType::AWAKE) << " " << moment->get_moment(MCType::ASLEEP) << "\n";
            };
        };
        
        // Close
        f.close();
    };
    
    void IxnParamTraj::write_adjoint_traj_to_file(int timepoint_start, int no_timesteps, std::string fname, bool with_timepoints) const {
        std::ofstream f;
        
        // Open
        f.open(fname);
        
        // Make sure we found it
        if (!f.is_open()) {
            std::cerr << ">>> Error: IxnParamTraj::write_moment_traj_to_file <<< could not write to file: " << fname << std::endl;
            exit(EXIT_FAILURE);
        };
        
        // Go through all time
        if (with_timepoints) {
            for (auto timepoint=timepoint_start; timepoint<=timepoint_start+no_timesteps; timepoint++) {
                f << timepoint << " " << _adjoint->get_val_at_timepoint(timepoint) << "\n";
            };
        } else {
            for (auto timepoint=timepoint_start; timepoint<=timepoint_start+no_timesteps; timepoint++) {
                f << _adjoint->get_val_at_timepoint(timepoint) << "\n";
            };
        };
        
        // Close
        f.close();
    };

};
