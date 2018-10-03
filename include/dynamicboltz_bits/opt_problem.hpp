#include <string>
#include <vector>

#ifndef FWDS_IXN_PARAM_H
#define FWDS_IXN_PARAM_H
#include "fwds/fwds_ixn_param.hpp"
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	// Forwards
	class IxnParam;
	class Lattice;

	/****************************************
	Filename timeseries
	****************************************/

	struct FNameSeries {

		std::vector<std::string> fnames;

	};

	/****************************************
	Filename collection
	****************************************/

	class FNameSeriesColl {

	private:

		// Filenames
		std::vector<FNameSeries> _fnames;
		std::vector<int> _idxs;

	public:

		/********************
		Get fnames
		********************/

		const std::vector<FNameSeries>& get_fname_series_all() const;
		const FNameSeries& get_fname_series(int idx) const;

		/********************
		Add fname
		********************/

		void add_fname_series(FNameSeries fname_series);
		void clear();

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
		bool VERBOSE_NU = false;
		bool VERBOSE_ADJOINT = false;
		bool VERBOSE_UPDATE = false;
		bool VERBOSE_WAKE_ASLEEP = false;
		bool VERBOSE_MOMENT = true;

		// Random integral
		bool MODE_random_integral_range = false;
		int VAL_random_integral_range_size = 10;

		// Nesterov
		bool nesterov = true;
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

		void init_structures(int no_timesteps, int batch_size);

		/********************
		Wake/asleep loop
		********************/

		// timepoint_start & timepoint_end = inclusive
		void wake_sleep_loop(int timepoint_start, int timepoint_end, int batch_size, int no_latt_sampling_steps, FNameSeriesColl &fname_coll, bool verbose=false);

		/********************
		Solve
		********************/

		// Check if options passed are valid
		void check_options(int no_timesteps, int batch_size, double dt, double dopt, int no_latt_sampling_steps, OptionsSolve options);

		// One step
		void solve_one_step(int i_opt_step, int no_timesteps, int batch_size, double dt, double dopt, int no_latt_sampling_steps, FNameSeriesColl &fname_coll, OptionsSolve options = OptionsSolve(), bool should_check_options=true);

		// Many steps
		void solve(int no_opt_steps, int no_timesteps, int batch_size, double dt, double dopt, int no_latt_sampling_steps, FNameSeriesColl &fname_coll, OptionsSolve options = OptionsSolve());


	};

};