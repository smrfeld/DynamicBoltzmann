#ifndef LIST_H
#define LIST_H
#include <list>
#endif

#ifndef STRING_H
#define STRING_H
#include <string>
#endif

#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

#define DIAG_SETUP 0
#define DIAG_SOLVE 0

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Struct for specifying a dimension
	****************************************/

	// Type of dimension
	enum DimType { H, J, W, B };

	struct Dim {
		// Name
		std::string name;

		// Type
		DimType type;

		// Name of associated species
		std::string species1;
		std::string species2;

		// Min/max/npts
		double min,max;
		int n;

		// Initial value
		double init;

		// Basis function dimension names
		std::vector<std::string> basis_func_dims;

		// Constructor
		Dim(std::string name, DimType type, std::string species, std::vector<std::string> basis_func_dims, double min, double max, int n, double init);
		Dim(std::string name, DimType type, std::string species1, std::string species2, std::vector<std::string> basis_func_dims, double min, double max, int n, double init);
	};

	/****************************************
	OptProblem
	****************************************/

	class OptProblem {

	private:

		class Impl;
		std::unique_ptr<Impl> _impl;

	public:

		/********************
		Constructor
		********************/

		/**
		 * @brief      Constructor
		 * @param[in]  dims         Vector of dimensions
		 * @param[in]  species      Vector of species present
		 * @param[in]  t_max        Maximum time
		 * @param[in]  n_t          No. of timepoints in the trajectory
		 * @param[in]  batch_size   Batch size
		 * @param[in]  box_length   Box length
		 * @param[in]  dopt         Optimization step interval
		 * @param[in]  n_opt        No. optimization steps
		 * @param[in]  lattice_dim  The lattice dimension
		 */
		OptProblem(std::vector<Dim> dims, std::vector<std::string> species, double t_max, int n_t, int batch_size, int box_length, double dopt, int n_opt, int lattice_dim=3);

		/**
		 * @brief      Move constructor (movable but no copies)
		 * @param[in]  other  The other
		 */
		OptProblem(OptProblem&& other);

	    /**
	     * @brief      Move assignment (movable but no copies)
	     * @param[in]  other  The other
	     * @return     Moved
	     */
	    OptProblem& operator=(OptProblem&& other);

	    /**
	     * @brief      Destructor
	     */
		~OptProblem();

		/********************
		Set properties	
		********************/

		// Any dim
		void add_hidden_unit(std::vector<std::vector<int>> lattice_idxs, std::string species);
		// 1D specific
		void add_hidden_unit(std::vector<int> lattice_idxs, std::string species);

		/**
		 * @brief      Sets the dir for i/o.
		 * @param[in]  dir   The dir
		 */
		void set_dir_io(std::string dir);

		/**
		 * @brief      Sets the filename index to start reading.
		 * @param[in]  idx   The index
		 */
		void set_fname_start_idx(int idx);
		
		/**
		 * @brief      Adds a filename for the lattices.
		 * @param[in]  f     The filename
		 */
		void add_fname(std::string f);

		// Set the number of CD steps
		void set_n_cd_steps(int n_steps);

		/********************
		Validate setup
		********************/

		/**
		 * @brief      Validate the setup by printing.
		 */
		void validate_setup() const;

		/********************
		Solve interaction parameter traj
		********************/

		void solve_ixn_param_traj();

		/********************
		Solve for variational trajectory
		********************/

		void solve_var_traj();

		/********************
		Solve
		********************/

		void solve(bool verbose=false, bool same_lattice=false);
		void solve_varying_ic(bool verbose=false);

		/********************
		Read some initial conditions
		********************/

		void read_init_cond(std::string dir);

		/********************
		Write
		********************/

		void write_bf_grids() const;
		void write_t_grid() const;

		void write_ixn_params(std::string dir, int idx) const;
		void write_ixn_params(std::string dir, int idx1, int idx2) const;
		void write_bfs(std::string dir, int idx) const;
		void write_var_terms(std::string dir, int idx) const;
		void write_moments(std::string dir, int idx) const;
		void write_moments(std::string dir, int idx1, int idx2) const;

		void set_flag_write_bf_only_final();

		/********************
		Read
		********************/

		void read_bf(std::string bf_name, std::string fname);
	};


};