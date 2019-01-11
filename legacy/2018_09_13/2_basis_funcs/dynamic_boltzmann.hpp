#ifndef BASISFUNC_h
#define BASISFUNC_h
#include "basis_func.hpp"
#endif

#ifndef MAP_h
#define MAP_h
#include <map>
#endif

#ifndef LIST_h
#define LIST_h
#include <list>
#endif

#ifndef SPECIES_h
#define SPECIES_h
#include "species.hpp"
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
		BasisFunc2D _f_h;
		BasisFunc2D _f_j;

		// The variational problem solution
		double* _h_traj;
		double* _j_traj;
		double*** _var_hh_traj;
		double*** _var_hj_traj;
		double*** _var_jh_traj;
		double*** _var_jj_traj;
		double* _t_grid;
		double* _h_grid;
		double* _j_grid;

		// Species present
		std::list<Species> _species;

		// Filenames to choose from
		std::vector<std::string> _fnames;

		// Nu range
		double _h_min,_h_max,_j_min,_j_max;

		// Number nu
		int _n_h,_n_j;

		// Increment
		double _dh,_dj;

		// Max time
		double _t_max;

		// Number timesteps
		int _n_t;

		// Timestep
		double _dt;

		// Number of steps in this nu solution
		int _n_t_soln;

		// Initial value for nu
		double _h_init,_j_init;

		// Delta function at some time index, nu value
		double _delta(double h1, double h2, double j1, double j2);

		// Batch size
		int _n_batch;

		// Number of annealing steps
		int _n_annealing;

		// Moments from data, annealing
		std::map<Species*,double> _moms_h_awake;
		std::map<Species*,double> _moms_h_asleep;
		std::map<Species*,std::map<Species*,double>> _moms_j_awake;
		std::map<Species*,std::map<Species*,double>> _moms_j_asleep;

		// Lattice size
		int _box_length;

		// Lattice to hold the current sample of the batch
		Lattice _latt;

		// Update step for optimization
		double _dopt;

		// Number opt steps
		int _n_opt;

	public:

		/********************
		Constructor
		********************/

		OptProblem(double h_min, double h_max, int n_h, double j_min, double j_max, int n_j, double t_max, int n_t, double h_init, double j_init, int batch_size, int n_annealing, int box_length, double dopt, int n_opt);
		~OptProblem();

		/********************
		Set properties	
		********************/

		void add_species(std::string sp);
		void add_fname(std::string f);

		/********************
		Solve for h,j
		********************/

		void solve_hj_traj();

		/********************
		Solve for variational trajectory
		********************/

		void solve_var_traj();

		/********************
		Write
		********************/

		void write_hj_traj(std::string fname);
		void write_var_traj(std::string fname);
		void write_bfs(std::string fname_h, std::string fname_j);

		/********************
		Solve
		********************/

		void solve(bool verbose=false);
	};

};