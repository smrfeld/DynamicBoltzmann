#ifndef LIST_h
#define LIST_h
#include <list>
#endif

#ifndef STRING_h
#define STRING_h
#include <string>
#endif

#ifndef VECTOR_h
#define VECTOR_h
#include <vector>
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

		Dim(std::string name, DimType type, double guess);
		Dim(std::string name, DimType type, std::string species, double guess);
		Dim(std::string name, DimType type, std::vector<std::string> species, double guess);
		Dim(std::string name, DimType type, std::vector<std::vector<std::string>> species, double guess);
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

		// Guess
		double guess() const;

		// Get species
		std::vector<std::string> get_species_h() const;
		std::vector<std::string> get_species_b() const;
		std::vector<std::vector<std::string>> get_species_J() const;
		std::vector<std::vector<std::string>> get_species_K() const;
		std::vector<std::vector<std::string>> get_species_W() const;

		/********************
		Setters
		********************/

		// Add species
		void add_species_h(std::string species);
		void add_species_b(std::string species);
		void add_species_J(std::string species1, std::string species2);
		void add_species_K(std::string species1, std::string species2, std::string species3);
		void add_species_W(std::string species_visible, std::string species_hidden);
	};

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
		Add counter
		********************/

		void add_counter(std::string s);
		void add_counter(std::string s1, std::string s2);
		void add_counter(std::string s1, std::string s2, std::string s3);
		void add_counter(std::string s1, std::string s2, std::string s3, std::string s4);
	};
};


