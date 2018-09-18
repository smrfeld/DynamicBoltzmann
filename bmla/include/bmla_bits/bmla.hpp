#ifndef STRING_H
#define STRING_H
#include <string>
#endif

#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

/************************************
* Namespace for bmla
************************************/

namespace bmla {

	// Forward
	class Dim;

	/****************************************
	Struct to specify choices for solving BMLA
	****************************************/

	struct OptionsSolveBMLA {

		// Are the following units binary/probabilistic
		bool awake_visible_are_binary;
		bool awake_hidden_are_binary;
		bool asleep_visible_are_binary;
		bool asleep_hidden_are_binary;
		bool asleep_final_visible_are_binary;
		bool asleep_final_hidden_are_binary;

		// Verbosity
		bool verbose;

		// If the MSE dips below this, quit
		bool mse_quit_mode;
		double mse_quit;

		// L2 Reg
		bool l2_reg_mode;
		double l2_lambda;

		// Limit l2 to certain parms
		std::vector<std::string> l2_limit_to_params;

		// Use a single lattice, irregardless of batch size
		bool use_single_lattice;

		// Filename to write the solution traj to
		bool write_soln_traj;
		std::string fname_write_soln_traj;

		// Filename to write the moment traj to
		bool write_moment_traj;
		std::string fname_write_moment_traj;

		// Track the solution traj
		bool track_soln_traj;

		// Nesterov mode
		bool nesterov;

		// Writing index
		int opt_idx_start_writing;

		// Whether to append to existing files
		bool append;

		// Whether to initialize from a filename, or from a random lattice
		bool start_CD_with_random_latt;

		/********************
		Constructor
		********************/

		OptionsSolveBMLA() {
			awake_visible_are_binary = true;
			awake_hidden_are_binary = true;
			asleep_visible_are_binary = true;
			asleep_hidden_are_binary = true;
			asleep_final_visible_are_binary = true;
			asleep_final_hidden_are_binary = true;
			verbose = true;
			mse_quit_mode = false;
			mse_quit = 0.0;
			l2_reg_mode = false;
			l2_lambda = 0.0;
			use_single_lattice = false;
			write_soln_traj = false;
			fname_write_soln_traj = "";
			write_moment_traj = false;
			fname_write_moment_traj = "";
			track_soln_traj = false;
			nesterov = true;
			opt_idx_start_writing = 0;
			append = false;
			start_CD_with_random_latt = false;
		};
	};

	/****************************************
	Struct to specify choices for sampling
	****************************************/

	struct OptionsSampling {

		// Are the following units binary/probabilistic
		bool awake_visible_are_binary;
		bool awake_hidden_are_binary;
		bool asleep_visible_are_binary;
		bool asleep_hidden_are_binary;
		bool asleep_final_visible_are_binary;
		bool asleep_final_hidden_are_binary;

		// Verbosity
		bool verbose;

		// What moments to report
		std::vector<std::string> report_counts;
		std::vector<std::pair<std::string,std::string>> report_nns;
		std::vector<std::vector<std::string>> report_triplets;
		std::vector<std::vector<std::string>> report_quartics;

		// Write
		// If batch_size = 1, writes the traj
		// If batch_size > 1, writes the average
		bool write_traj;
		std::string fname_write_traj;

		/********************
		Constructor
		********************/

		OptionsSampling() {
			awake_visible_are_binary = true;
			awake_hidden_are_binary = true;
			asleep_visible_are_binary = true;
			asleep_hidden_are_binary = true;
			asleep_final_visible_are_binary = true;
			asleep_final_hidden_are_binary = true;
			verbose = true;
			write_traj = false;
			fname_write_traj = "";
		};
	};

	/****************************************
	BMLA
	****************************************/

	class BMLA {

	private:

		class Impl;
		std::unique_ptr<Impl> _impl;

	public:

		/********************
		Constructor
		********************/

		BMLA(std::vector<Dim> dims, std::vector<std::string> species_visible, std::vector<std::string> species_hidden, int box_length, int lattice_dim=3);
		BMLA(const BMLA& other);
		BMLA(BMLA&& other);
		BMLA& operator=(BMLA other);
		~BMLA();

		/********************
		Set a parameter value for dim
		********************/

		void set_param_for_dim(std::string dim_name, double val);

		/********************
		Set a fixed awake moment
		********************/

		void set_fixed_awake_moment_for_dim(std::string dim_name, double val);

		/********************
		Add a hidden unit
		********************/

		// Any dim
		void add_hidden_unit(std::vector<std::string> species_possible, std::vector<std::vector<int>> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params);
		// 1D specific
		void add_hidden_unit(std::vector<std::string> species_possible, std::vector<int> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params);

		// Validate
		void validate_hidden() const;

		/********************
		Solve for the h,j corresponding to a given lattice
		********************/

		void solve(std::string fname, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolveBMLA options=OptionsSolveBMLA());
		void solve(std::vector<std::string> fnames, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolveBMLA options=OptionsSolveBMLA());

		/********************
		At the current ixns params, sample and report the specified moments
		********************/

		void sample(int batch_size, int n_cd_steps, OptionsSampling options=OptionsSampling());

		/********************
		Read param file
		********************/

		void read(std::string fname);

		/********************
		Write
		********************/

		void write(std::string fname, bool append=false);
		void write(std::string fname, int idx, bool append=false);
		void write_ave(std::string fname, int last_n_steps, bool append=false);
		void write_ave(std::string fname, int last_n_steps, int idx, bool append=false);

		/********************
		Write lattice
		********************/

		void write_lattice(std::string fname);

	};
};

