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
    class FNameColl;

	/****************************************
	OptProblemRBM Options
	****************************************/

	enum class Solver : unsigned int { SGD, NESTEROV, ADAM };

	struct OptionsWakeSleepRBM {

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

		// Parallel
		bool parallel = true;
	};

	struct OptionsSolveRBM {
		// Should check options before starting
		bool should_check_options = true;

		// Verbosity
		bool verbose_update = false;
		bool verbose_moment = true;

		// L2 Reg mode
		bool l2_reg = false;
		std::map<std::shared_ptr<IxnParam>,double> l2_lambda;
		std::map<std::shared_ptr<IxnParam>,double> l2_center;

		// Variable learning rate
		bool var_learning_rates = false;
		std::map<std::shared_ptr<IxnParam>,double> var_learning_rate_values;

		// Options for Wake-Sleep loop
		OptionsWakeSleepRBM options_wake_sleep = OptionsWakeSleepRBM();

		// Options for the solvers
		Solver solver = Solver::ADAM;

		// Nesterov
		double nesterov_acc = 0.5;

		// Adam
		double adam_beta_1 = 0.9;
		double adam_beta_2 = 0.999;
		double adam_eps = 0.00000001;
	};

	/****************************************
	OptProblemRBM
	****************************************/

	class OptProblemRBM {

	private:

		// Ixn params
		std::vector<std::shared_ptr<IxnParam>> _ixn_params;

		// Lattice
		std::shared_ptr<Lattice> _latt;

		// Constructor helpers
		void _clean_up();
		void _move(OptProblemRBM &other);
		void _copy(const OptProblemRBM& other);

	public:

		/********************
		Constructor
		********************/

		OptProblemRBM(std::shared_ptr<Lattice> latt, std::vector<std::shared_ptr<IxnParam>> ixn_params);
		OptProblemRBM(const OptProblemRBM& other);
		OptProblemRBM(OptProblemRBM&& other);
		OptProblemRBM& operator=(const OptProblemRBM &other);
		OptProblemRBM& operator=(OptProblemRBM &&other);
		~OptProblemRBM();

		/********************
		Init structures
		********************/

		void init_structures(int batch_size);

		/********************
		Wake/asleep loop
		********************/

		void wake_sleep_loop(int batch_size, int no_cd_sampling_steps, FNameColl &fname_coll, OptionsWakeSleepRBM options=OptionsWakeSleepRBM());

		/********************
		Solve
		********************/

		// Check if options passed are valid
		void check_options(int batch_size, double dopt, int no_cd_sampling_steps, OptionsSolveRBM options);

		// One step
		void solve_one_step(int i_opt_step, int batch_size, double dopt, int no_cd_sampling_steps, FNameColl &fname_coll, OptionsSolveRBM options = OptionsSolveRBM());

		// Many steps
		void solve(int no_opt_steps, int batch_size, double dopt, int no_cd_sampling_steps, FNameColl &fname_coll, OptionsSolveRBM options = OptionsSolveRBM());
	};

};
