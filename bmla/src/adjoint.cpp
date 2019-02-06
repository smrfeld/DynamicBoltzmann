#include "../include/dblz_bits/adjoint.hpp"

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

	/****************************************
	Adjoint
	****************************************/

    // ***************
    // MARK: - Constructor
    // ***************

	Adjoint::Adjoint(std::string name, ITptr ixn_param_traj) {
		_name = name;
		_ixn_param_traj = ixn_param_traj;

        set_no_timesteps(0);
	};
	Adjoint::Adjoint(const Adjoint& other) {
		_copy(other);
	};
	Adjoint::Adjoint(Adjoint&& other) {
		_move(other);
	};
    Adjoint& Adjoint::operator=(const Adjoint& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    Adjoint& Adjoint::operator=(Adjoint&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };
	Adjoint::~Adjoint()
	{
		_clean_up();
	};
	void Adjoint::_clean_up() {
        // Nothing....
    };
	void Adjoint::_copy(const Adjoint& other) {
		_name = other._name;
		_ixn_param_traj = other._ixn_param_traj;
		_no_timesteps = other._no_timesteps;
		_no_timepoints = other._no_timepoints;
        _vals = other._vals;
		_timepoint_zero_end_cond = other._timepoint_zero_end_cond;
	};
	void Adjoint::_move(Adjoint& other) {
		_name = other._name;
        _ixn_param_traj = other._ixn_param_traj;
		_no_timesteps = other._no_timesteps;
		_no_timepoints = other._no_timepoints;
        _vals = other._vals;
		_timepoint_zero_end_cond = other._timepoint_zero_end_cond;

		// Reset the other
		other._name = "";
		other._ixn_param_traj = nullptr;
		other._no_timesteps = 0;
		other._no_timepoints = 0;
        other._vals.clear();
        other._timepoint_zero_end_cond = 0;
	};

    // ***************
    // MARK: - Timesteps
    // ***************
    
	int Adjoint::get_no_timesteps() const {
		return _no_timesteps;
	};
	void Adjoint::set_no_timesteps(int no_timesteps) {
		_no_timesteps = no_timesteps;
		_no_timepoints = _no_timesteps + 1;
        
        while (_vals.size() < _no_timepoints) {
            _vals.push_back(0.0);
        };
        while (_vals.size() > _no_timepoints) {
            _vals.pop_back();
        };
	};

    // ***************
    // MARK: - Init cond
    // ***************
    
	int Adjoint::get_timepoint_zero_end_cond() const {
		return _timepoint_zero_end_cond;
	};
	void Adjoint::set_timepoint_zero_end_cond(int timepoint) {
		if (timepoint >= _no_timepoints) {
			std::cout << ">>> Error: Adjoint::set_timepoint_zero_end_cond <<< no timepoints is: " << _no_timepoints << std::endl;
			exit(EXIT_FAILURE);
		};

		_timepoint_zero_end_cond = timepoint;
		_vals[_timepoint_zero_end_cond] = 0.0;
	};
    
    // ***************
    // MARK: - Name, type
    // ***************

	std::string Adjoint::get_name() const {
		return _name;
	};

	ITptr Adjoint::get_ixn_param_traj() const {
		return _ixn_param_traj;
	};

    // ***************
    // MARK: - Value
    // ***************
    
	double Adjoint::get_val_at_timepoint(int timepoint) const {
		if (timepoint >= _no_timepoints) {
			std::cerr << ">>> Error: Adjoint::get_val_at_timepoint <<< " << timepoint << " is out of bounds: " << _no_timepoints << std::endl;
			exit(EXIT_FAILURE);
		};
		return _vals[timepoint];
	};

    // ***************
    // MARK: - Diff eq
    // ***************
    
	void Adjoint::solve_diff_eq_at_timepoint_to_minus_one(int timepoint, double dt, bool l2_mode, const std::map<ITptr,double> &l2_lambda, const std::map<ITptr,double> &l2_center) {
		if (timepoint >= _no_timepoints) {
			std::cerr << ">>> Error: Adjoint::solve_diff_eq_at_timepoint_to_minus_one <<< " << timepoint << " is out of bounds: " << _no_timepoints << std::endl;
			exit(EXIT_FAILURE);
		};
		if (timepoint > _timepoint_zero_end_cond) {
			std::cerr << ">>> Error: Adjoint::solve_diff_eq_at_timepoint_to_minus_one <<< " << timepoint << " is beyond the zero endpoint: " << _timepoint_zero_end_cond << std::endl;
			exit(EXIT_FAILURE);		
		};

		// Deriv val
		double deriv_term=0.0;
		double deriv, adjoint_val;

		// Go through dependencies (numerator F and adjoint)
        for (auto const dep_pair: _ixn_param_traj->get_diff_eq_dependencies()) {

			// Derivative
			deriv = dep_pair.first->get_deriv_wrt_nu_at_timepoint(timepoint,dep_pair.second);

			// Adjoint
			auto adjoint = dep_pair.first->get_parent_ixn_param_traj()->get_adjoint();
			if (!adjoint) {
				std::cerr << ">>> Error: Adjoint::solve_diff_eq_at_timepoint_to_minus_one <<< No adjoint for ixn param: " << dep_pair.first->get_parent_ixn_param_traj()->get_name() << std::endl;
				exit(EXIT_FAILURE);	
			};
			adjoint_val = adjoint->get_val_at_timepoint(timepoint);

			// Add
			deriv_term += deriv * adjoint_val;

		};

		// Get domain
		/*
		auto diff_eq_rhs = _ixn_param->get_diff_eq_rhs();
		if (!diff_eq_rhs) {
			std::cerr << ">>> Error: Adjoint::solve_diff_eq_at_timepoint_to_minus_one <<< No diff eq rhs for ixn param!" << std::endl;
			exit(EXIT_FAILURE);	
		};

		int i_dim = 0;
		for (auto const &domain: diff_eq_rhs->get_domain()) {
			// Get ixn param
			auto ixn_param_in_deriv = domain->get_ixn_param();

			// Get adjoint
			auto adjoint_in_deriv = ixn_param_in_deriv->get_adjoint();
			if (!adjoint_in_deriv) {
				std::cerr << ">>> Error: Adjoint::solve_diff_eq_at_timepoint_to_minus_one <<< No adjoint for ixn param: " << ixn_param_in_deriv->get_name() << std::endl;
				exit(EXIT_FAILURE);	
			};

			// Val
			deriv += diff_eq_rhs->get_deriv_wrt_nu_at_timepoint(timepoint,i_dim) * adjoint_in_deriv->get_val_at_timepoint(timepoint);

			// Next dim
			i_dim++;
		}; 
		*/

		// Difference in moments
        auto moment = _ixn_param_traj->get_ixn_param_at_timepoint(timepoint)->get_moment();
        double moment_delta = moment->get_moment(MCType::ASLEEP) - moment->get_moment(MCType::AWAKE);
        
		// L2 reg
		if (l2_mode) {
			// L2 mode

            double ixn_param_val = _ixn_param_traj->get_ixn_param_at_timepoint(timepoint)->get_val();
			double l2_term = 2.0 * l2_lambda.at(_ixn_param_traj) * (ixn_param_val-l2_center.at(_ixn_param_traj));
			std::cout << _name << " " << timepoint << " " << moment_delta << " " << l2_term << std::endl;

			// Step
			_vals[timepoint-1] = _vals[timepoint] - dt * (moment_delta - deriv + l2_term);

		} else {
			// Not L2 mode

			// Step
			_vals[timepoint-1] = _vals[timepoint] - dt * (moment_delta - deriv);
		};
	};

	/********************
	Write to file
	********************/

	void Adjoint::write_to_file(std::string fname) const {
		std::ofstream f;

		// Open
		f.open(fname);

		// Make sure we found it
		if (!f.is_open()) {
			std::cerr << ">>> Error: Adjoint::write_to_file <<< could not write to file: " << fname << std::endl;
			exit(EXIT_FAILURE);
		};

		// Go through all time
		for (auto timepoint=0; timepoint<_no_timepoints; timepoint++) {
			f << _vals[timepoint] << "\n";
		};

		// Close
		f.close();
	};
};
