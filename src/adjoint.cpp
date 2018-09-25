#include "../include/dynamicboltz_bits/adjoint.hpp"

// Other headers
#include "../include/dynamicboltz_bits/ixn_param.hpp"
#include "../include/dynamicboltz_bits/general.hpp"
#include "../include/dynamicboltz_bits/diff_eq_rhs.hpp"
#include "../include/dynamicboltz_bits/moment.hpp"

#include <iostream>
#include <fstream>

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	Adjoint - IMPLEMENTATION DEFINITIONS
	****************************************/

	/********************
	Constructor
	********************/

	Adjoint::Adjoint(std::string name, Iptr ixn_param) {
		_name = name;
		_ixn_param = ixn_param;

		_no_timesteps = 0;
		_no_timepoints = 1;
		_zero_end_cond_timepoint = 0;

		_vals = new double[_no_timepoints];
		std::fill_n(_vals,_no_timepoints,0.0);
		_vals[0] = 0.0;
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
		safeDelArr(_vals);
	};
	void Adjoint::_copy(const Adjoint& other) {
		_name = other._name;
		_ixn_param = other._ixn_param;
		_no_timesteps = other._no_timesteps;
		_no_timepoints = other._no_timepoints;
		_vals = new double[_no_timepoints];
		std::copy(other._vals,other._vals+_no_timepoints,_vals);
		_zero_end_cond_timepoint = other._zero_end_cond_timepoint;
	};
	void Adjoint::_move(Adjoint& other) {
		_name = other._name;
		_ixn_param = other._ixn_param;
		_no_timesteps = other._no_timesteps;
		_no_timepoints = other._no_timepoints;
		_vals = new double[_no_timepoints];
		std::copy(other._vals,other._vals+_no_timepoints,_vals);
		_zero_end_cond_timepoint = other._zero_end_cond_timepoint;

		// Reset the other
		other._name = "";
		other._ixn_param = nullptr;
		other._no_timesteps = 0;
		other._no_timepoints = 0;
		safeDelArr(other._vals);
		other._zero_end_cond_timepoint = 0;
	};

	/********************
	Validate setup
	********************/

	void Adjoint::check_setup() const {
		// ...
	};

	/********************
	Timesteps
	********************/

	int Adjoint::get_no_timesteps() const {
		return _no_timesteps;
	};
	void Adjoint::set_no_timesteps(int no_timesteps) {
		_no_timesteps = no_timesteps;
		_no_timepoints = _no_timesteps + 1;
		// Vals
		safeDelArr(_vals);
		_vals = new double[_no_timepoints];
		std::fill_n(_vals,_no_timepoints,0.0);
		_vals[_no_timepoints-1] = 0.0;
	};


	/********************
	Init cond
	********************/

	int Adjoint::get_zero_end_cond_timepoint() const {
		return _zero_end_cond_timepoint;
	};
	void Adjoint::set_zero_end_cond_timepoint(int timepoint) {
		if (timepoint >= _no_timepoints) {
			std::cout << ">>> Error: Adjoint::set_zero_end_cond_timepoint <<< no timepoints is: " << _no_timepoints << std::endl;
			exit(EXIT_FAILURE);
		};

		_zero_end_cond_timepoint = timepoint;
		_vals[_zero_end_cond_timepoint] = 0.0;
	};

	/********************
	Name, type
	********************/

	std::string Adjoint::get_name() const {
		return _name;
	};

	Iptr Adjoint::get_ixn_param() const {	
		return _ixn_param;
	};

	/********************
	Value
	********************/

	double Adjoint::get_val_at_timepoint(int timepoint) const {
		if (timepoint >= _no_timepoints) {
			std::cerr << ">>> Error: Adjoint::get_val_at_timepoint <<< " << timepoint << " is out of bounds: " << _no_timepoints << std::endl;
			exit(EXIT_FAILURE);
		};
		return _vals[timepoint];
	};

	/********************
	Diff eq
	********************/

	void Adjoint::solve_diff_eq_at_timepoint_to_minus_one(int timepoint, double dt) {
		if (timepoint >= _no_timepoints) {
			std::cerr << ">>> Error: Adjoint::solve_diff_eq_at_timepoint_to_minus_one <<< " << timepoint << " is out of bounds: " << _no_timepoints << std::endl;
			exit(EXIT_FAILURE);
		};
		if (timepoint > _zero_end_cond_timepoint) {
			std::cerr << ">>> Error: Adjoint::solve_diff_eq_at_timepoint_to_minus_one <<< " << timepoint << " is beyond the zero endpoint: " << _zero_end_cond_timepoint << std::endl;
			exit(EXIT_FAILURE);		
		};

		auto diff_eq_rhs = _ixn_param->get_diff_eq_rhs();
		if (!diff_eq_rhs) {
			std::cerr << ">>> Error: Adjoint::solve_diff_eq_at_timepoint_to_minus_one <<< No diff eq rhs for ixn param!" << std::endl;
			exit(EXIT_FAILURE);	
		};

		// Deriv val
		double deriv=0.0;

		// Get domain
		int i_dim = 0;
		for (auto const &domain: diff_eq_rhs->get_domain()) {
			// Get ixn param
			auto ixn_param_in_deriv = domain.get_ixn_param();

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

		// Difference in moments
		double moment_delta = _ixn_param->get_moment()->get_moment_at_timepoint(MomentType::ASLEEP, timepoint) - _ixn_param->get_moment()->get_moment_at_timepoint(MomentType::AWAKE, timepoint);

		// Step
		_vals[timepoint-1] = _vals[timepoint] - dt * (deriv - moment_delta);
	};

	/********************
	Reset to zero
	********************/

	void Adjoint::reset_to_zero() {
		for (auto timepoint=0; timepoint<_no_timepoints; timepoint++) {
			_vals[timepoint] = 0.0;
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