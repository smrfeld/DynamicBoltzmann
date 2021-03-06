#include "../include/dblz_bits/adjoint.hpp"

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
        
        _timepoint_zero_end_cond = _no_timepoints-1;
        _vals[_timepoint_zero_end_cond] = 0.0;
	};

    // ***************
    // MARK: - Init cond
    // ***************
    
    /*
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
    */
    
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
    // MARK: - Set zero endpoint
    // ***************
    
    void Adjoint::set_timepoint_zero_end_cond(int timepoint) {
        _timepoint_zero_end_cond = timepoint;
        _vals[_timepoint_zero_end_cond] = 0.0;
    };

    // ***************
    // MARK: - Value
    // ***************
    
	double Adjoint::get_val_at_timepoint(int timepoint) const {
		if (timepoint >= _no_timepoints) {
            std::cerr << ">>> Error: Adjoint::get_val_at_timepoint <<< " << timepoint << " is out of bounds: no timesteps: " << _no_timesteps << " no timepoints: " << _no_timepoints << std::endl;
			exit(EXIT_FAILURE);
		};
		return _vals[timepoint];
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
