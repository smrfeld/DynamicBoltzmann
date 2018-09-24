#ifndef STRING_H
#define STRING_H
#include <string>
#endif

#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
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

		/********************
		Get random subset
		********************/

		std::vector<int> get_random_subset(int size);
	};

	/****************************************
	OptProblem Options
	****************************************/

	struct OptionsSolve {
		bool VERBOSE_NU = false;
		bool VERBOSE_ADJOINT = false;
		bool VERBOSE_UPDATE = false;
		bool VERBOSE_WAKE_ASLEEP = false;
		bool VERBOSE_MOMENT = true;
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

		void wake_asleep_loop(int no_timesteps, int batch_size, int no_latt_sampling_steps, FNameSeriesColl &fname_coll, bool verbose=false);

		/********************
		Solve
		********************/

		// One step
		void solve_one_step(int no_timesteps, int batch_size, double dt, double dopt, int no_latt_sampling_steps, FNameSeriesColl &fname_coll, OptionsSolve options = OptionsSolve());

		// Many steps
		void solve(int no_opt_steps, int no_timesteps, int batch_size, double dt, double dopt, int no_latt_sampling_steps, FNameSeriesColl &fname_coll, OptionsSolve options = OptionsSolve());


	};

};