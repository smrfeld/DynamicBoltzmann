#ifndef STRING_H
#define STRING_H
#include <string>
#endif

#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

#ifndef MAP_H
#define MAP_H
#include <map>
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	Filename structure - only needed for solve_varying_ic
	****************************************/

	struct FName {

		// Filename
		std::string fname;

		// Idxs, possibly more than one
		int idx_split;
		int idx_sample;

		// Initial condition
		std::string fname_ic;

		// Whether to write this idx
		bool write;

		// Constructor
		FName(std::string fname, int idx_split, int idx_sample, std::string fname_ic, bool write=true) {
			this->fname=fname;
			this->idx_split=idx_split;
			this->idx_sample=idx_sample;
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

		// Write ixn params
		bool write_ixn_params;

		// Write moments
		bool write_moments;

		// Write only the final basis function
		bool write_bf_only_final;

		// L2 Reg for parameters
		bool l2_reg_params_mode;
		std::map<std::string,double> l2_lambda_params;
		std::map<std::string,double> l2_reg_centers;

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

		// Exponential decay for the time
		// exp( - lamdbda * ( t - t0 ) )
		// where t=0....# timesteps (not real time!)
		bool exp_decay;
		std::vector<double> exp_decay_t0_values;
		std::vector<double> exp_decay_lambda_values;

		// Time cutoffs
		bool time_cutoff;
		std::vector<int> time_cutoff_start_values; // inclusive
		std::vector<int> time_cutoff_end_values; // exclusive

		/********************
		Constructor
		********************/

		OptionsSolve() {
			verbose = true;
			nesterov = true;
			use_same_lattice_in_batch = false;
			write = false;
			dir_write = "";
			write_ixn_params = true;
			write_moments = true;
			write_bf_only_final = false;
			l2_reg_params_mode = false;
			clear_dir = true;
			time_idx_start_reading = 0;
			opt_idx_start_writing = 0;
			awake_visible_are_binary = true;
			awake_hidden_are_binary = true;
			asleep_visible_are_binary = true;
			asleep_hidden_are_binary = true;
			asleep_final_visible_are_binary = true;
			asleep_final_hidden_are_binary = true;
			exp_decay = false;
			time_cutoff = false;
		};
	};

	/****************************************
	OptProblem
	****************************************/

	// Forward
	class Dim;

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
		 * @param[in]  no_timepoints          No. of timepoints in the trajectory
		 * @param[in]  batch_size   Batch size
		 * @param[in]  box_length   Box length
		 * @param[in]  dopt         Optimization step interval
		 * @param[in]  n_opt        No. optimization steps
		 * @param[in]  lattice_dim  The lattice dimension
		 */
		OptProblem(std::vector<Dim> dims, std::vector<std::string> species_visible, std::vector<std::string> species_hidden, double t_max, int no_timepoints, int box_length, int lattice_dim=3);
		OptProblem(const OptProblem& other);
		OptProblem(OptProblem&& other);
	    OptProblem& operator=(OptProblem other);
		~OptProblem();

		/********************
		Change the time limit
		********************/

		void set_no_timepoints(int no_timepoints);

		/********************
		Set IC for ixn param
		********************/

		void set_ic_for_ixn_param(std::string param_name, double val);

		/********************
		Set a fixed awake moment
		********************/

		void set_fixed_awake_moment_for_dim(std::string param_name, std::vector<double> vals);

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
		void validate_graph() const;

		/********************
		Solve interaction parameter traj
		********************/

		void solve_ixn_param_traj();

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

		// Single IC, dividing up randomly the traj and restarting the diff eq integration for the var terms
		void solve_rand_div(std::vector<std::string> fnames, int n_opt, int batch_size, int n_cd_steps, double dopt, int n_divisions, OptionsSolve options=OptionsSolve());
		void solve_rand_div(std::vector<std::vector<std::string>> fnames, int n_opt, int n_cd_steps, double dopt, int n_divisions, OptionsSolve options=OptionsSolve());

		/********************
		Read basis function
		********************/

		void read_basis_func(std::string bf_name, std::string fname);

		void read_init_cond_for_ixn_param(std::string ixn_func_name, std::string fname, int line_idx=-1);

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
		void write_moments(std::string dir, int idx) const;
		void write_moments(std::string dir, int idx1, int idx2) const;
		void write_moments(std::string dir, int idx, std::vector<int> idxs) const;

	};

};