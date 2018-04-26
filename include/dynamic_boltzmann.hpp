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

#ifndef UTILITY_H
#define UTILITY_H
#include <utility>
#endif

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
	Filename structure - only needed for solve_varying_ic
	****************************************/

	struct FName {

		// Filename
		std::string fname;

		// Idxs, possibly more than one
		std::vector<int> idxs;

		// Initial condition
		std::string fname_ic;

		// Whether to write this idx
		bool write;

		FName(std::string fname, int idx, std::string fname_ic, bool write=true) : FName(fname,std::vector<int>({idx}),fname_ic,write) {};
		FName(std::string fname, std::vector<int> idxs, std::string fname_ic, bool write=true) {
			this->fname=fname;
			this->idxs=idxs;
			this->fname_ic=fname_ic;
			this->write = write;
		};
	};

	/****************************************
	Options for solving
	****************************************/

	struct OptionsSolve {

		// Verbose
		bool verbose;

		// Nesterov
		bool nesterov;

		// Use the same lattice for all samples in the batch
		bool use_same_lattice_in_batch;

		// Whether to write data
		bool write;
		std::string dir_write;

		// Write the variational terms
		bool write_var_terms;

		// Write only the final basis function
		bool write_bf_only_final;

		// Clear directory
		bool clear_dir;

		// Time indexes to start
		int time_idx_start_reading;

		// Optimization step to start writing
		int opt_idx_start_writing;

		// Are the following units binary/probabilistic
		bool awake_visible_are_binary;
		bool awake_hidden_are_binary;
		bool asleep_visible_are_binary;
		bool asleep_hidden_are_binary;
		bool asleep_final_visible_are_binary;
		bool asleep_final_hidden_are_binary;

		// Locality factor for variational term
		bool local_decay;
		double local_decay_factor;

		/********************
		Constructor
		********************/

		OptionsSolve() {
			verbose = true;
			nesterov = true;
			use_same_lattice_in_batch = false;
			write = false;
			dir_write = "";
			write_var_terms = false;
			write_bf_only_final = false;
			clear_dir = true;
			time_idx_start_reading = 0;
			opt_idx_start_writing = 0;
			awake_visible_are_binary = true;
			awake_hidden_are_binary = true;
			asleep_visible_are_binary = true;
			asleep_hidden_are_binary = true;
			asleep_final_visible_are_binary = true;
			asleep_final_hidden_are_binary = true;
			local_decay = false;
			local_decay_factor = 0.;
		};
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
		OptProblem(std::vector<Dim> dims, double t_max, int n_t, int box_length, int lattice_dim=3);

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
		Set hidden layer	
		********************/

		// Any dim
		void add_hidden_unit(std::vector<std::vector<int>> lattice_idxs, std::string species);
		// 1D specific
		void add_hidden_unit(std::vector<int> lattice_idxs, std::string species);
		
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

		void solve(std::vector<std::string> fnames, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolve options=OptionsSolve());
		void solve_varying_ic(std::vector<FName> fnames, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolve options=OptionsSolve());
		void solve_varying_ic(std::vector<FName> fnames, std::vector<FName> fnames_used_in_every_batch, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolve options=OptionsSolve());

		/********************
		Read basis function
		********************/

		void read_basis_func(std::string bf_name, std::string fname);

		/********************
		Write
		********************/

		// Grids
		void write_bf_grids(std::string dir) const;
		void write_t_grid(std::string dir) const;

		// Solved values
		void write_ixn_params(std::string dir, int idx) const;
		void write_ixn_params(std::string dir, int idx1, int idx2) const;
		void write_ixn_params(std::string dir, int idx, std::vector<int> idxs) const;		
		void write_bfs(std::string dir, int idx) const;
		void write_var_terms(std::string dir, int idx) const;
		void write_moments(std::string dir, int idx) const;
		void write_moments(std::string dir, int idx1, int idx2) const;
		void write_moments(std::string dir, int idx, std::vector<int> idxs) const;

	};

};