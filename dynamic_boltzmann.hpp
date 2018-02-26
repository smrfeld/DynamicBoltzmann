#ifndef BASISFUNC_h
#define BASISFUNC_h
#include "basis_func_nd.hpp"
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

// Diagnostic flags
#define DIAG_INIT 0
#define DIAG_VAR 0
#define DIAG_SOLVE 0

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	OptProblem
	****************************************/

	class OptProblem {

	private:

		// Number of dimensions
		int _n_dim;

		// Vector of dimensions present
		std::vector<Dim> _dims;

		// The basis function solutions
		// Size: (1+_n_dim)
		//    Dim 1 = basis funcs, 1 to _n_dim
		//    Dim 2 to _n_dim+1 = the grid
		// (Vector b/c there is no default constructor)
		std::vector<BasisFuncND> _bfs;

		// Grids to hold the solutions
		// Size: 2
		//    Dim 1 = time
		//    Dim 2 = values, 1 to _n_dim
		double **_soln_traj;

		// Grids to hold the variational problem solutions
		//  Size: (3+_n_dim)
		//     Dim 1 = time
		//     Dim 2 = numerator = delta nu to vary at (_n_dim possible)
		//     Dim 3 = denominator = delta f to vary with respect to (_n_dim possible)
		//     Dim 4 to _n_dim+3 = grid points for the denominator
		// (Vector b/c there is no default constructor)
		std::vector<std::vector<std::vector<GridND>>> _var_traj;

		// Updates to the basis functions
		// Size: (1+_n_dim)
		//    Dim 1 = basis funcs, 1 to _n_dim
		//    Dim 2 to _n_dim+1 = the grid
		// (Vector b/c there is no default constructor)
		std::vector<GridND> _dbfs;

		// Initial points
		std::vector<double> _init;

		// Time dimension
		Dim _time;

		// Species present
		std::list<Species> _species;

		// Filenames to choose from
		std::vector<std::string> _fnames;

		// Number of steps in this nu solution
		int _n_t_soln;

		// The current time in the optimization step
		int _t_opt;

		// Batch size
		int _n_batch;

		// Number of annealing steps
		int _n_annealing;

		// Lattice size
		int _box_length;

		// Lattice to hold the current sample of the batch
		Lattice _latt;

		// Update step for optimization
		double _dopt;

		// Number opt steps
		int _n_opt;

		// Recursive function to iterate over all grid points for the variational term
		void _iterate_var(int dim);
		// Indexes for the recursion function
		// Size: 1
		//     Dim 1: grid values, 1 to _n_dim
		int *_var_idxs;
		// Time, num, denom idx
		int _var_time, _var_nu, _var_f;
		// Derivatives of the basis functions at the solution corresponding to the current timepoint
		// Size: 2
		//     Dim 1: basis functions, 1 to _n_dim
		//     Dim 2: diff. variable, 1 to _n_dim
		double **_bfs_derivs_var;
		// Delta function
		double _var_delta();

		// Find a species by name
		Species* _find_species_by_name(std::string name);
		// Find a interaction param by Species
		int _find_param_by_species(Species *sp);
		int _find_param_by_species(Species *sp1, Species *sp2);

	public:

		/********************
		Constructor
		********************/

		OptProblem(int n_dim, double t_max, int n_t, int batch_size, int n_annealing, int box_length, double dopt, int n_opt, std::vector<Dim> dims, std::vector<double> init, std::vector<std::string> species);
		~OptProblem();

		/********************
		Set properties	
		********************/

		void add_fname(std::string f);

		/********************
		Solve for solution traj
		********************/

		void solve_traj();

		/********************
		Solve for variational trajectory
		********************/

		void solve_var_traj();

		/********************
		Write
		********************/

		void write_soln_traj(std::string fname);
		void write_var_traj(std::string fname);
		void write_bfs(std::string fname);
		void write_bfs(std::string dir, int idx);
		void write_grid(std::string fname);
		void write_moments(std::string fname, bool append);

		/********************
		Solve
		********************/

		void solve(bool verbose=false);
	};

};