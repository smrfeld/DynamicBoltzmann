#include "dynamic_boltzmann.hpp"
#include <numeric> 
#include "general.hpp"

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	OptProblem
	****************************************/

	/********************
	Constructor
	********************/

	OptProblem::OptProblem(double nu_min, double nu_max, int n_nu, double t_max, int n_t, double nu_init) : _f(nu_min,nu_max,n_nu)
	{
		this->_t_max = t_max;
		this->_n_t = n_t;
		this->_n_t_soln = n_t;
		this->_dt = t_max / (n_t - 1);
		this->_nu_init = nu_init;

		// All nu_init by default
		std::iota(this->_nu_traj.begin(), this->_nu_traj.end(), nu_init);

		// Time values
		for (int i=0; i<n_t; i++) {
			this->_tvals.push_back(i*this->_dt);
		};

		// TMP
		_f.test_fill();
	};
	OptProblem::~OptProblem() {};

	/********************
	Solve for nu traj
	********************/

	void OptProblem::solve_nu_traj() 
	{
		_nu_traj.clear();
		_tvals.clear();
		_nu_traj.push_back(_nu_init);
		_tvals.push_back(0.);
		double n;
		for (int i=0; i<_n_t-1; i++) {
			n = _nu_traj[i] + _dt * _f.get_val(_nu_traj[i]);
			if (n <= _f.get_nu_max() && n >= _f.get_nu_min()) {
				_nu_traj.push_back(n);
				_tvals.push_back((i+1)*_dt);
			} else {
				_n_t_soln = i+1;
				return;
			};
		};
		_n_t_soln = _n_t;
	};

	/********************
	Print
	********************/
	void OptProblem::print_nu_traj() {
		for (int i=0; i<_n_t_soln; i++) 
		{ 
			std::cout << _tvals[i] << " " << _nu_traj[i] << std::endl;
		};
	};

};

