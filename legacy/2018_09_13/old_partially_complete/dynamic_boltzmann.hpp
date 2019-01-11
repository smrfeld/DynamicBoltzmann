
#ifndef LATTICE_h
#define LATTICE_h
#include "lattice.hpp"
#endif

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	OptProblem
	****************************************/

	class OptProblem {

	private:

		// The solution
		BasisFunc _f;

		// The variational problem solution
		std::vector<double> _nu_traj;
		std::vector<Lattice> _var_traj;
		std::vector<double> _t_grid;
		Lattice _nu_grid;

		// Max time
		double _t_max;

		// Number timesteps
		int _n_t;

		// Timestep
		double _dt;

		// Number of steps in this nu solution
		int _n_t_soln;

		// Initial value for nu
		double _nu_init;

		// Delta function at some time index, nu value
		double _delta(double nu1, double nu2);

	public:

		// Constructor
		OptProblem(double nu_min, double nu_max, int n_nu, double t_max, int n_t, double nu_init);
		// Destructor
		~OptProblem();

		// Solve for nu
		void solve_nu_traj();

		// Solve for variational trajectory
		void solve_var_traj();

		// Print nu solution
		void print_nu_traj();
		void print_var_traj();

		// Write the variational traj

	};

};