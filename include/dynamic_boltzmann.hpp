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
	enum DimType { H, J, K, W, B };

	class Dim {

	private:

		class Impl;
		std::unique_ptr<Impl> _impl;

	public:

		/********************
		Constructor
		********************/

		Dim(std::string name, DimType type, std::vector<std::string> basis_func_dims, double min, double max, int n, double init);
		Dim(std::string name, DimType type, std::string species, std::vector<std::string> basis_func_dims, double min, double max, int n, double init);
		Dim(std::string name, DimType type, std::vector<std::string> species, std::vector<std::string> basis_func_dims, double min, double max, int n, double init);
		Dim(std::string name, DimType type, std::vector<std::vector<std::string>> species, std::vector<std::string> basis_func_dims, double min, double max, int n, double init);
		Dim(const Dim& other);
		Dim(Dim&& other);
		Dim& operator=(Dim other);
		~Dim();

		/********************
		Getters
		********************/

		// Name
		std::string name() const;

		// Type
		DimType type() const;

		// Does it apply to any species?
		bool any_species() const;

		// Basis func dims
		std::vector<std::string> basis_func_dims() const;

		// Min/max/n/init
		double min() const;
		double max() const;
		double n() const;
		double init() const;

		// Get species
		std::vector<std::string> get_species_h() const;
		std::vector<std::string> get_species_b() const;
		std::vector<std::vector<std::string>> get_species_J() const;
		std::vector<std::vector<std::string>> get_species_K() const;
		std::vector<std::vector<std::string>> get_species_W() const;

		/********************
		Setters
		********************/

		// Add basis func dimension
		void add_basis_func_dim(std::string dim);

		// Add species
		void add_species_h(std::string species);
		void add_species_b(std::string species);
		void add_species_J(std::string species1, std::string species2);
		void add_species_K(std::string species1, std::string species2, std::string species3);
		void add_species_W(std::string species_visible, std::string species_hidden);
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

		// Constructor
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
		OptProblem(std::vector<Dim> dims, std::vector<std::string> species_visible, std::vector<std::string> species_hidden, double t_max, int n_t, int box_length, int lattice_dim=3);
		OptProblem(const OptProblem& other);
		OptProblem(OptProblem&& other);
	    OptProblem& operator=(OptProblem other);
		~OptProblem();

		/********************
		Set hidden layer	
		********************/

		// Any dim
		void add_hidden_unit(std::vector<std::string> species_possible, std::vector<std::vector<int>> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params);
		// 1D specific
		void add_hidden_unit(std::vector<std::string> species_possible, std::vector<int> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params);
		
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

		// Single IC
		void solve(std::vector<std::string> fnames, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolve options=OptionsSolve());
		void solve(std::vector<std::vector<std::string>> fname_collection, int n_opt, int n_cd_steps, double dopt, OptionsSolve options=OptionsSolve());

		// Varying IC
		void solve_varying_ic(std::vector<FName> fnames, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolve options=OptionsSolve());
		void solve_varying_ic(std::vector<FName> fnames, std::vector<FName> fnames_used_in_every_batch, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolve options=OptionsSolve());
		void solve_varying_ic(std::vector<std::vector<FName>> fname_collection, int n_opt, int n_cd_steps, double dopt, OptionsSolve options=OptionsSolve());

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