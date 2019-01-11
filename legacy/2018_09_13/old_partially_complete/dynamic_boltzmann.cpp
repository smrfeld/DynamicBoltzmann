#include "dynamic_boltzmann.hpp"
#include <numeric> 
#include "general.hpp"
#include "math.h"

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

	OptProblem::OptProblem(double nu_min, double nu_max, int n_nu, double t_max, int n_t, double nu_init) : _f(nu_min,nu_max,n_nu), _nu_grid(nu_min,nu_max,n_nu)
	{
		this->_t_max = t_max;
		this->_n_t = n_t;
		this->_n_t_soln = n_t;
		this->_dt = t_max / (n_t - 1);
		this->_nu_init = nu_init;

		// Fill the _nu_grid incrementally
	    for (int i=0; i<n_nu; i++) {
	    	_nu_grid.vals[i] = nu_min + i*_nu_grid.delta;
	    };

		// TEMP: fill the basis func with sin
	    for (int i=0; i<n_t; i++) {
			_f.vals[i] = sin(2.0*M_PI*i/(n_t-1))+0.1;
		};
		_f.update_derivs();
	};
	OptProblem::~OptProblem() {};

	/********************
	Solve for nu traj
	********************/

	void OptProblem::solve_nu_traj() 
	{
		_nu_traj.clear();
		_t_grid.clear();
		_nu_traj.push_back(_nu_init);
		_t_grid.push_back(0.);
		double n;
		for (int i=0; i<_n_t-1; i++) {
			n = _nu_traj[i] + _dt * _f.get_val(_nu_traj[i]);
			if (n <= _f.max && n >= _f.min) {
				_nu_traj.push_back(n);
				_t_grid.push_back((i+1)*_dt);
			} else {
				_n_t_soln = i+1;
				return;
			};
		};
		_n_t_soln = _n_t;
	};

	/********************
	Solve for var traj
	********************/

	void OptProblem::solve_var_traj()
	{
		// Clear the old
		_var_traj.clear();

		// Initial
		_var_traj.push_back(Lattice(_f.min,_f.max,_f.n));
		// (all zero by default)

		// Go through all times
		for (int i=0; i<_n_t_soln-1; i++) {
			_var_traj.push_back(Lattice(_f.min,_f.max,_f.n));
			std::cout << "Old" << std::endl;
			std::cout << _var_traj[i] << std::endl;
			// Go through all nu
			for (int j=0; j<_f.n; j++) {
				//std::cout << _nu_grid[j] << " " << _nu_traj[i] << " delta = " << _delta(_nu_traj[i],_nu_grid[j]) << std::endl;
				_var_traj[i+1][j] = _var_traj[i][j] + _dt * ( _delta(_nu_traj[i],_nu_grid[j]) + _f.get_deriv(_nu_traj[i]) * _var_traj[i][j] );
				// std::cout << _var_traj[i+1][j] << std::endl;
			};
			std::cout << "New" << std::endl;
			std::cout << _var_traj[i+1] << std::endl;
		};
	};

	/********************
	Print
	********************/

	void OptProblem::print_nu_traj() {
		for (int i=0; i<_n_t_soln; i++) 
		{ 
			std::cout << _t_grid[i] << " " << _nu_traj[i] << std::endl;
		};
	};
	void OptProblem::print_var_traj() {
		for (int i=0; i<_n_t_soln; i++) 
		{
			for (int j=0; j<_f.n; j++) {
				std::cout << _t_grid[i] << " " << _nu_grid[j] << " " << _var_traj[i][j] << std::endl;
			};
		};
	};


	/****************************************
	OptProblem - PRIVATE
	****************************************/

	/********************
	Delta function at some time index, nu value
	********************/

	double OptProblem::_delta(double nu1, double nu2) 
	{
		return exp(-pow(nu1-nu2,2)/(2.*0.01)) / sqrt(2. * M_PI * 0.01);
	};



};

