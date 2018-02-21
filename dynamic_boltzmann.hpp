#ifndef BASISFUNC_h
#define BASISFUNC_h
#include "basis_func.hpp"
#endif

#ifndef MAP_h
#define MAP_h
#include <map>
#endif

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
		double* _nu_traj;
		double** _var_traj;
		double* _t_grid;
		double* _nu_grid;

		// Species present
		std::vector<std::string> _species;

		// Nu range
		double _nu_min,_nu_max;

		// Number nu
		int _n_nu;

		// Increment
		double _dnu;

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

		// Batch size
		int _n_batch;

		// Number of annealing steps
		int _n_annealing;

		// Moments from data, annealing
		std::map<std::string,double> _moms_awake;
		std::map<std::string,double> _moms_asleep;

		// Lattice size
		int _box_length;

		// Data batch
		Lattice* _data;

		// Update step for optimization
		double _dopt;

	public:

		// Constructor
		OptProblem(double nu_min, double nu_max, int n_nu, double t_max, int n_t, double nu_init, int batch_size, int n_annealing, std::vector<std::string> species_list, int box_length, double dopt);
		// Destructor
		~OptProblem();

		// Solve for nu
		void solve_nu_traj();

		// Solve for variational trajectory
		void solve_var_traj();

		// Print nu solution
		void print_nu_traj();
		void print_var_traj();

		// Write
		void write_nu_traj(std::string fname);
		void write_var_traj(std::string fname);
		void write_bfs(std::string fname);

		// Run
		void solve();
	};

};