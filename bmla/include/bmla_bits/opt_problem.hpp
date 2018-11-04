#include <string>
#include <vector>
#include <memory>

/************************************
* Namespace for bmla
************************************/

namespace bmla {

	// Forwards
	class IxnParam;
	class Lattice;

	/****************************************
	Filename collection
	****************************************/

	class FNameColl {

	private:

		// Filenames
		std::vector<std::string> _fnames;
		std::vector<int> _idxs;

	public:

		/********************
		Get fnames
		********************/

		const std::vector<std::string>& get_fnames_all() const;
		const std::string& get_fname(int idx) const;

		/********************
		Add fname
		********************/

		void add_fname(std::string fname);

		/********************
		Get random subset
		********************/

		std::vector<int> get_random_subset(int size);
	};

	/****************************************
	OptProblem Options
	****************************************/

	struct OptionsSolve {
		// Verbosity
		bool VERBOSE_UPDATE = false;
		bool VERBOSE_WAKE_ASLEEP = false;
		bool VERBOSE_MOMENT = true;

		// Start with random lattice
		bool start_with_random_lattice = false;
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

		void wake_sleep_loop(int batch_size, int no_latt_sampling_steps, FNameColl &fname_coll, bool verbose=false, bool start_with_random_lattice=false);

		/********************
		Solve
		********************/

		// Check if options passed are valid
		void check_options(int batch_size, double dopt, int no_latt_sampling_steps, OptionsSolve options);

		// One step
		void solve_one_step(int batch_size, double dopt, int no_latt_sampling_steps, FNameColl &fname_coll, OptionsSolve options = OptionsSolve(), bool should_check_options=true);

		// Many steps
		void solve(int no_opt_steps, int batch_size, double dopt, int no_latt_sampling_steps, FNameColl &fname_coll, OptionsSolve options = OptionsSolve());
	};

};
