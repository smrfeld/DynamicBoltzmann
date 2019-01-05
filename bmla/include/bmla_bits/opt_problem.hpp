#include <string>
#include <vector>
#include <memory>
#include <map>

/************************************
* Namespace for bmla
************************************/

namespace bmla {

	// Forwards
	class IxnParam;
	class Lattice;

	/****************************************
	Filename
	****************************************/

	struct FName {

		std::string name;
		bool binary;

		FName(std::string name, bool binary);
	};

	/****************************************
	Filename collection
	****************************************/

	class FNameColl {

	private:

		// Filenames
		std::vector<FName> _fnames;
		std::vector<int> _idxs;

	public:

		/********************
		Get fnames
		********************/

		const std::vector<FName>& get_fnames_all() const;
		const FName& get_fname(int idx) const;

		/********************
		Add fname
		********************/

		void add_fname(FName fname);

		/********************
		Get random subset
		********************/

		std::vector<int> get_random_subset(int size);
	};

	/****************************************
	OptProblem Options
	****************************************/

	struct OptionsWakeSleep {

		// Verbosity
		bool verbose = false;

		// Start with random lattice
		bool start_with_random_lattice = false;

		// Awake phase
		// Is the hidden layer binary for the awake moment?
		bool is_awake_moment_hidden_binary = false;
		// Binarize hidden layer after awake moment?
		// (if !awake_moment_hidden_binary)
		bool should_binarize_hidden_after_awake_moment = true;

		// Asleep phase
		// Is the visible reconstruction binary?
		bool is_asleep_visible_binary = true;
		// Is the hidden reconstruction binary, EXCEPT in the last phase?
		bool is_asleep_hidden_binary = true;
		// Is the hidden reconstruction binary in the last phase?
		bool is_asleep_hidden_binary_final = false;

		// Layerwise sampling
		bool layer_wise = true;

		// Parallel
		bool parallel = true;
	};

	struct OptionsSolve {
		// Should check options before starting
		bool should_check_options = true;

		// Verbosity
		bool VERBOSE_UPDATE = false;
		bool VERBOSE_MOMENT = true;

		// L2 Reg mode
		bool MODE_l2_reg = false;
		std::map<std::shared_ptr<IxnParam>,double> VAL_l2_lambda;
		std::map<std::shared_ptr<IxnParam>,double> VAL_l2_center;

		// Variable learning rate
		bool MODE_var_learning_rates = false;
		std::map<std::shared_ptr<IxnParam>,double> VAL_var_learning_rates;

		// Nesterov
		bool nesterov = true;

		// Options for Wake-Sleep loop
		OptionsWakeSleep options_wake_sleep = OptionsWakeSleep();
	};

	/****************************************
	OptProblem
	****************************************/

	class OptProblem {

	private:

		// Ixn params
		std::vector<std::shared_ptr<IxnParam>> _ixn_params;

		// Lattice
		std::shared_ptr<Lattice> _latt;

		// Constructor helpers
		void _clean_up();
		void _move(OptProblem &other);
		void _copy(const OptProblem& other);

	public:

		/********************
		Constructor
		********************/

		OptProblem(std::shared_ptr<Lattice> latt);
		OptProblem(const OptProblem& other);
		OptProblem(OptProblem&& other);
		OptProblem& operator=(const OptProblem &other);
		OptProblem& operator=(OptProblem &&other);
		~OptProblem();

		/********************
		Setters
		********************/

		void add_ixn_param(std::shared_ptr<IxnParam> ixn_param);
		void set_lattice(std::shared_ptr<Lattice> latt);

		/********************
		Init structures
		********************/

		void init_structures(int batch_size);

		/********************
		Wake/asleep loop
		********************/

		void wake_sleep_loop(int batch_size, int no_latt_sampling_steps, FNameColl &fname_coll, OptionsWakeSleep options=OptionsWakeSleep());

		/********************
		Solve
		********************/

		// Check if options passed are valid
		void check_options(int batch_size, double dopt, int no_latt_sampling_steps, OptionsSolve options);

		// One step
		void solve_one_step(int i_opt_step, int batch_size, double dopt, int no_latt_sampling_steps, FNameColl &fname_coll, OptionsSolve options = OptionsSolve());

		// Many steps
		void solve(int no_opt_steps, int batch_size, double dopt, int no_latt_sampling_steps, FNameColl &fname_coll, OptionsSolve options = OptionsSolve());
	};

};
