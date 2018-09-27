#include "../include/dynamicboltz_bits/opt_problem.hpp"

// Other headers
#include "../include/dynamicboltz_bits/ixn_param.hpp"
#include "../include/dynamicboltz_bits/lattice.hpp"
#include "../include/dynamicboltz_bits/adjoint.hpp"
#include "../include/dynamicboltz_bits/moment.hpp"
#include "../include/dynamicboltz_bits/diff_eq_rhs.hpp"
#include "../include/dynamicboltz_bits/general.hpp"

#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	Filename collection
	****************************************/

	/********************
	Get fnames
	********************/

	const std::vector<FNameSeries>& FNameSeriesColl::get_fname_series_all() const {
		return _fnames;
	};
	const FNameSeries& FNameSeriesColl::get_fname_series(int idx) const {
		return _fnames[idx];
	};

	/********************
	Add fname
	********************/

	void FNameSeriesColl::add_fname_series(FNameSeries fname_series) {
		_fnames.push_back(fname_series);
		_idxs.push_back(_fnames.size()-1);
	};

	/********************
	Get random subset
	********************/

	std::vector<int> FNameSeriesColl::get_random_subset(int size) {
		if (_fnames.size() < size) {
			std::cerr << ">>> Error: FNameSeriesColl::get_random_subset <<< size of subset is equal to size of filename collection" << std::endl;
			exit(EXIT_FAILURE);
		};

		// Shuffle idxs
	    std::random_device rd;
	    std::mt19937 g(rd());
	    std::shuffle(_idxs.begin(), _idxs.end(), g);

	    std::vector<int>::const_iterator first = _idxs.begin();
	    std::vector<int>::const_iterator last = _idxs.begin() + size;
	    return std::vector<int>(first,last);
	};

	/****************************************
	OptProblem
	****************************************/
	
	/********************
	Constructor
	********************/

	OptProblem::OptProblem(std::shared_ptr<Lattice> latt) {
		_latt = latt;
	};
	OptProblem::OptProblem(const OptProblem& other) {
		_copy(other);
	};
	OptProblem::OptProblem(OptProblem&& other) {
		_move(other);
	};
	OptProblem& OptProblem::operator=(const OptProblem& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	OptProblem& OptProblem::operator=(OptProblem&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
	};
	OptProblem::~OptProblem() {
		_clean_up();
	};

	void OptProblem::_clean_up() {};
	void OptProblem::_move(OptProblem &other) {
		_ixn_params = std::move(other._ixn_params);
		_latt = std::move(other._latt);
	};
	void OptProblem::_copy(const OptProblem& other) {
		_ixn_params = other._ixn_params;
		_latt = other._latt;
	};

	/********************
	Setters
	********************/

	void OptProblem::add_ixn_param(std::shared_ptr<IxnParam> ixn_param) {
		_ixn_params.push_back(ixn_param);
	};
	void OptProblem::set_lattice(std::shared_ptr<Lattice> latt) {
		_latt = latt;
	};

	/********************
	Init structures
	********************/

	void OptProblem::init_structures(int no_timesteps, int batch_size) {

		// Ixn params
		for (auto &ixn_param: _ixn_params) {
			ixn_param->set_no_timesteps(no_timesteps);

			// Adjoint
			if (ixn_param->get_adjoint()) {
				ixn_param->get_adjoint()->set_no_timesteps(no_timesteps);
				ixn_param->get_adjoint()->set_zero_end_cond_timepoint(no_timesteps);
			};

			// Moment
			ixn_param->get_moment()->set_no_timesteps(no_timesteps);
			ixn_param->get_moment()->set_batch_size(batch_size);
		};
	};

	/********************
	Wake/asleep loop
	********************/

	void OptProblem::wake_sleep_loop(int timepoint_start, int timepoint_end, int batch_size, int no_latt_sampling_steps, FNameSeriesColl &fname_coll, bool verbose) {
		if (verbose) {
			std::cout << "--- Sampling lattice ---" << std::endl;
		};

		// Make a subset
		std::vector<int> idx_subset = fname_coll.get_random_subset(batch_size);

		// Reset all moments
		for (auto &ixn_param: _ixn_params) {
			ixn_param->get_moment()->reset_to_zero();
		};

		// Iterate over time
		for (int timepoint=timepoint_start; timepoint<=timepoint_end; timepoint++) {

			if (verbose) {
				std::cout << "(" << timepoint << " in " << timepoint_start << " ... " << timepoint_end << ") : " << std::flush;
			};

			// Iterate over the batch
			for (int i_batch=0; i_batch<batch_size; i_batch++)
			{

				if (verbose) {
					std::cout << "." << std::flush;
				};

				// Convert lattice to binary mode
				_latt->all_units_convert_to_b_mode();

				// Read latt
				_latt->read_from_file(fname_coll.get_fname_series(idx_subset[i_batch]).fnames[timepoint]); // binary units

				// Sample hidden
				_latt->sample_h(); // binary units

				// Reap awake
				for (auto &ixn_param: _ixn_params) {
					ixn_param->get_moment()->reap_as_timepoint_in_batch(MomentType::AWAKE, timepoint, i_batch);
				};

				// Convert lattice to prob mode
				_latt->all_units_convert_to_p_mode();

				// Sample
				for (int i_sampling_step=0; i_sampling_step<no_latt_sampling_steps; i_sampling_step++) 
				{
					_latt->sample_v_at_timepoint(timepoint,false); // prob units
					_latt->sample_h_at_timepoint(timepoint,false); // prob units
				};

				// Print
				// latt.print_occupancy();

				// Reap asleep
				for (auto &ixn_param: _ixn_params) {
					ixn_param->get_moment()->reap_as_timepoint_in_batch(MomentType::ASLEEP, timepoint, i_batch, false);
				};
			};
			
			// Average moments at this timepoint
			for (auto &ixn_param: _ixn_params) {
				ixn_param->get_moment()->average_reaps_as_timepoint(MomentType::AWAKE, timepoint);
				ixn_param->get_moment()->average_reaps_as_timepoint(MomentType::ASLEEP, timepoint);
			};

			if (verbose) {
				std::cout << std::endl;
			};

		};

		if (verbose) {
			std::cout << "--- [Finished] Sampled lattice ---" << std::endl;
			std::cout << std::endl;
		};
	};

	/********************
	Solve
	********************/

	// Check if options passed are valid
	void OptProblem::check_options(int no_timesteps, int batch_size, double dt, double dopt, int no_latt_sampling_steps, OptionsSolve options) {

		if (options.MODE_random_integral_range) {
			if (options.VAL_random_integral_range_size > no_timesteps) {
				std::cerr << ">>> Error: OptProblem::check_options <<< options must satisfy: VAL_random_integral_range_size <= no timesteps!" << std::endl;
				exit(EXIT_FAILURE);
			};
		};

	};

	// One step
	void OptProblem::solve_one_step(int no_timesteps, int batch_size, double dt, double dopt, int no_latt_sampling_steps, FNameSeriesColl &fname_coll, OptionsSolve options, bool should_check_options) {

		int no_timepoints = no_timesteps + 1;

		/*****
		Check options
		*****/

		if (should_check_options) {
			check_options(no_timesteps,batch_size,dt,dopt,no_latt_sampling_steps,options);
		};

		/*****
		Option: random integrand range mode
		*****/

		int timepoint_integral_start=0;
		int timepoint_integral_end=no_timesteps;

		if (options.MODE_random_integral_range) {
			timepoint_integral_start = randI(0,no_timesteps-options.VAL_random_integral_range_size);
			timepoint_integral_end = timepoint_integral_start+options.VAL_random_integral_range_size;
		};

		/*****
		Solve diff eq for F
		*****/

		if (options.VERBOSE_NU) {
			std::cout << "--- Solving diff eq ---" << std::endl;
		};

		// Solve over all time
		for (auto t=0; t<no_timesteps; t++) {
			for (auto &ixn_param: _ixn_params) {
				ixn_param->solve_diff_eq_at_timepoint_to_plus_one(t,dt);
			};
		};

		if (options.VERBOSE_NU) {
			std::cout << "--- [Finished] Solving diff eq ---" << std::endl;
			std::cout << std::endl;
		};

		/*****
		Wake/asleep loop
		*****/

		wake_sleep_loop(timepoint_integral_start,timepoint_integral_end,batch_size,no_latt_sampling_steps,fname_coll,options.VERBOSE_WAKE_ASLEEP);

		if (options.VERBOSE_MOMENT) {
			for (auto &ixn_param: _ixn_params) {
				std::cout << ixn_param->get_name() << " " << std::flush;
				ixn_param->get_moment()->print_moment_comparison();
			};
		};

		/********************
		Solve diff eq for adjoint
		********************/

		if (options.VERBOSE_ADJOINT) {
			std::cout << "--- Solving adjoint ---" << std::endl;
		};

		// Reset
		for (auto &ixn_param: _ixn_params) {
			ixn_param->get_adjoint()->reset_to_zero();
		};

		for (auto t=timepoint_integral_end; t>=timepoint_integral_start+1; t--) {
			for (auto &ixn_param: _ixn_params) {
				ixn_param->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one(t,dt);
			};
		};

		if (options.VERBOSE_ADJOINT) {
			std::cout << "--- [Finished] Solving adjoint ---" << std::endl;
			std::cout << std::endl;
		};

		/********************
		Form the update
		********************/

		if (options.VERBOSE_UPDATE) {
			std::cout << "--- Calculating update ---" << std::endl;
		};

		for (auto &ixn_param: _ixn_params) {
			ixn_param->get_diff_eq_rhs()->update_calculate_and_store(timepoint_integral_start,timepoint_integral_end,dt,dopt);
		};

		if (options.VERBOSE_UPDATE) {
			std::cout << "--- [Finished] Calculating update ---" << std::endl;
			std::cout << std::endl;
		};

		/********************
		Committ the update
		********************/

		if (options.VERBOSE_UPDATE) {
			std::cout << "--- Committing update ---" << std::endl;
		};

		for (auto &ixn_param: _ixn_params) {
			ixn_param->get_diff_eq_rhs()->update_committ_stored();
		};

		if (options.VERBOSE_UPDATE) {
			std::cout << "--- [Finished] Committing update ---" << std::endl;
			std::cout << std::endl;
		};
	};

	// Many steps
	void OptProblem::solve(int no_opt_steps, int no_timesteps, int batch_size, double dt, double dopt, int no_latt_sampling_steps, FNameSeriesColl &fname_coll, OptionsSolve options) {

		/*****
		Init structures
		*****/

		init_structures(no_timesteps,batch_size);

		/*****
		Go...
		*****/

		for (int i_opt_step=1; i_opt_step<=no_opt_steps; i_opt_step++)
		{

			std::cout << "------------------" << std::endl;
			std::cout << "Opt step: " << i_opt_step << " / " << no_opt_steps << std::endl;
			std::cout << "------------------" << std::endl;

			// Solve
			solve_one_step(no_timesteps,batch_size,dt,dopt,no_latt_sampling_steps,fname_coll,options);
		};

	};

};