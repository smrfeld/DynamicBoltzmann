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

// Timing only
//#include <cstdio>
//#include <ctime>

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	Filename
	****************************************/

	FName::FName(std::string name, bool binary) {
		this->name = name;
		this->binary = binary;
	};

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
	void FNameSeriesColl::clear() {
		_fnames.clear();
		_idxs.clear();
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
		_moments.push_back(ixn_param->get_moment());
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

	void OptProblem::wake_sleep_loop(int timepoint_start, int timepoint_end, int batch_size, int no_latt_sampling_steps, FNameSeriesColl &fname_coll, OptionsWakeSleep options) {
		if (options.verbose) {
			std::cout << "--- Sampling lattice ---" << std::endl;
		};

		// Make a subset
		std::vector<int> batch_idx_subset = fname_coll.get_random_subset(batch_size);

		// Reset all moments
		for (auto &moment: _moments) {
			moment->reset_to_zero();
		};

		// Iterate over time
		for (int timepoint=timepoint_start; timepoint<=timepoint_end; timepoint++) {

			if (options.verbose) {
				std::cout << "(" << timepoint << " in " << timepoint_start << " ... " << timepoint_end << ") : " << std::flush;
			};

			// Iterate over the batch
			for (int i_batch=0; i_batch<batch_size; i_batch++)
			{

				if (options.verbose) {
					std::cout << "." << std::flush;
				};

				// clock_t t3 = clock();    

				// Read latt
				auto file = fname_coll.get_fname_series(batch_idx_subset[i_batch])[timepoint];
				_latt->read_from_file(file.name, file.binary);

				// clock_t t4 = clock();    

				// Sample hidden
				// Hidden: prob (= false)
				_latt->sample_up_v_to_h_at_timepoint(timepoint,options.layer_wise,options.is_awake_moment_hidden_binary);

				// clock_t t5 = clock();    

				// Reap awake
				for (auto &moment: _moments) {
					moment->reap_as_timepoint_in_batch(MomentType::AWAKE, timepoint, i_batch);
				};

				// clock_t t6 = clock();    

				// Convert hidden units to be binary
				if (!options.is_awake_moment_hidden_binary && options.should_binarize_hidden_after_awake_moment) {
					_latt->all_units_h_binarize();
				};

				// clock_t t7 = clock();    

				// Sample
				for (int i_sampling_step=0; i_sampling_step<no_latt_sampling_steps; i_sampling_step++) 
				{
					// Sample down (hidden -> visible)
					// Visible: binary (= true)
					// Hiddens: binary (= true)
					_latt->sample_down_h_to_v_at_timepoint(timepoint,options.layer_wise,options.is_asleep_visible_binary,options.is_asleep_hidden_binary);

					// Sample up (visible -> hidden)
					// If not last step:
					// Hiddens: binary (= true)
					// Else:
					// Hiddens: prob (= false)
					if (i_sampling_step != no_latt_sampling_steps-1) {
						_latt->sample_up_v_to_h_at_timepoint(timepoint,options.layer_wise,options.is_asleep_hidden_binary); // binary hiddens
					} else {
						_latt->sample_up_v_to_h_at_timepoint(timepoint,options.layer_wise,options.is_asleep_hidden_binary_final); // prob hiddens
					};
				};

				// clock_t t8 = clock();    

				// Print
				// latt.print_occupancy();

				// Reap asleep
				for (auto &moment: _moments) {
					moment->reap_as_timepoint_in_batch(MomentType::ASLEEP, timepoint, i_batch);
				};

				// clock_t t9 = clock();    

				// Timing
				/*
				double int1 = ( t4 - t3 ) / (double) CLOCKS_PER_SEC;
				double int2 = ( t5 - t4 ) / (double) CLOCKS_PER_SEC;
				double int3 = ( t6 - t5 ) / (double) CLOCKS_PER_SEC;
				double int4 = ( t7 - t6 ) / (double) CLOCKS_PER_SEC;
				double int5 = ( t8 - t7 ) / (double) CLOCKS_PER_SEC;
				double int6 = ( t9 - t8 ) / (double) CLOCKS_PER_SEC;
				double int_tot = int1 + int2 + int3 + int4 + int5 + int6;
				*/
				//std::cout << "Timepoint: " << timepoint << " batch: " << i_batch << std::endl;
				//std::cout << "Timing: " << int1 << " " << int2 << " " << int3 << " " << int4 << " " << int5 << " " << int6 << std::endl;
				//std::cout << "Percentages: " << int1/int_tot << " " << int2/int_tot << " " << int3/int_tot << " " << int4/int_tot << " " << int5/int_tot << " " << int6/int_tot << std::endl;
			};
			
			// Average moments at this timepoint
			for (auto &moment: _moments) {
				moment->average_reaps_as_timepoint(MomentType::AWAKE, timepoint);
				moment->average_reaps_as_timepoint(MomentType::ASLEEP, timepoint);
			};

			if (options.verbose) {
				std::cout << std::endl;
			};

		};

		if (options.verbose) {
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
			if (options.VAL_random_integral_range_end > no_timesteps) {
				std::cerr << ">>> Error: OptProblem::check_options <<< options must satisfy: VAL_random_integral_range_end <= no_timesteps!" << std::endl;
				exit(EXIT_FAILURE);
			};
		};
	};

	// One step
	void OptProblem::solve_one_step(int i_opt_step, int no_timesteps, int batch_size, double dt, double dopt, int no_latt_sampling_steps, FNameSeriesColl &fname_coll, OptionsSolve options) {

		int no_timepoints = no_timesteps + 1;

		/*****
		Check options
		*****/

		if (options.should_check_options) {
			check_options(no_timesteps,batch_size,dt,dopt,no_latt_sampling_steps,options);
		};

		/*****
		Option: random integrand range mode
		*****/

		// clock_t t1 = clock();    

		int timepoint_integral_start=0;
		int timepoint_integral_end=no_timesteps;

		if (options.MODE_random_integral_range) {
			timepoint_integral_start = randI(options.VAL_random_integral_range_start,options.VAL_random_integral_range_end);
			timepoint_integral_end = timepoint_integral_start+options.VAL_random_integral_range_size;
		};

		// clock_t t2 = clock();    
		// std::cout << "options: " << ( t2 - t1 ) / (double) CLOCKS_PER_SEC << std::endl;

		/*****
		Solve diff eq for F
		*****/

		if (options.VERBOSE_NU) {
			std::cout << "--- Solving diff eq ---" << std::endl;
		};

		// Solve over all time
		for (auto t=0; t<no_timesteps; t++) {
			for (auto &ixn_param: _ixn_params) {
				if (!ixn_param->get_is_val_fixed_to_init_cond() && !ixn_param->get_are_vals_fixed()) {
					ixn_param->solve_diff_eq_at_timepoint_to_plus_one(t,dt);
				};
			};
		};

		if (options.VERBOSE_NU) {
			std::cout << "--- [Finished] Solving diff eq ---" << std::endl;
			std::cout << std::endl;
		};

		// clock_t t3 = clock();    
		// std::cout << "diff eq F " << ( t3 - t2 ) / (double) CLOCKS_PER_SEC << std::endl;

		/*****
		Wake/asleep loop
		*****/

		clock_t t3 = clock();    

		wake_sleep_loop(timepoint_integral_start,timepoint_integral_end,batch_size,no_latt_sampling_steps,fname_coll,options.options_wake_sleep);

		if (options.VERBOSE_MOMENT) {
			for (auto &moment: _moments) {
				std::cout << moment->get_name() << " " << std::flush;
				moment->print_moment_comparison();
			};
		};

		clock_t t4 = clock();    
		std::cout << "wake sleep " << ( t4 - t3 ) / (double) CLOCKS_PER_SEC << std::endl;

		/********************
		Solve diff eq for adjoint
		********************/

		if (options.VERBOSE_ADJOINT) {
			std::cout << "--- Solving adjoint ---" << std::endl;
		};

		// Reset
		for (auto &ixn_param: _ixn_params) {
			if (!ixn_param->get_is_val_fixed_to_init_cond() && !ixn_param->get_are_vals_fixed()) {
				ixn_param->get_adjoint()->reset_to_zero();
			};
		};

		if (options.MODE_l2_reg) {
			for (auto t=timepoint_integral_end; t>=timepoint_integral_start+1; t--) {
				for (auto &ixn_param: _ixn_params) {
					if (!ixn_param->get_is_val_fixed_to_init_cond() && !ixn_param->get_are_vals_fixed()) {
						ixn_param->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one(t,dt,true,options.VAL_l2_lambda,options.VAL_l2_center);
					};
				};
			};
		} else {
			for (auto t=timepoint_integral_end; t>=timepoint_integral_start+1; t--) {
				for (auto &ixn_param: _ixn_params) {
					if (!ixn_param->get_is_val_fixed_to_init_cond() && !ixn_param->get_are_vals_fixed()) {
						ixn_param->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one(t,dt);
					};
				};
			};
		};

		if (options.VERBOSE_ADJOINT) {
			std::cout << "--- [Finished] Solving adjoint ---" << std::endl;
			std::cout << std::endl;
		};

		// clock_t t5 = clock();    
		// std::cout << "adjoint " << ( t5 - t4 ) / (double) CLOCKS_PER_SEC << std::endl;

		/********************
		Form the update
		********************/

		if (options.VERBOSE_UPDATE) {
			std::cout << "--- Calculating update ---" << std::endl;
		};

		double dopt_use = dopt;
		for (auto &ixn_param: _ixn_params) {
			if (!ixn_param->get_is_val_fixed_to_init_cond() && !ixn_param->get_are_vals_fixed()) {
				// Learning rate
				if (options.MODE_var_learning_rates) {
					dopt_use = options.VAL_var_learning_rates[ixn_param];
				};

				// Update
				ixn_param->get_diff_eq_rhs()->update_calculate_and_store(timepoint_integral_start,timepoint_integral_end,dt,dopt_use);
			};
		};

		if (options.VERBOSE_UPDATE) {
			std::cout << "--- [Finished] Calculating update ---" << std::endl;
			std::cout << std::endl;
		};

		// clock_t t6 = clock();    
		// std::cout << "form update " << ( t6 - t5 ) / (double) CLOCKS_PER_SEC << std::endl;

		/********************
		Committ the update
		********************/

		if (options.VERBOSE_UPDATE) {
			std::cout << "--- Committing update ---" << std::endl;
		};

		for (auto &ixn_param: _ixn_params) {
			if (!ixn_param->get_is_val_fixed_to_init_cond() && !ixn_param->get_are_vals_fixed()) {
				ixn_param->get_diff_eq_rhs()->update_committ_stored(options.nesterov,options.nesterov_acc);
			};
		};

		if (options.VERBOSE_UPDATE) {
			std::cout << "--- [Finished] Committing update ---" << std::endl;
			std::cout << std::endl;
		};

		// clock_t t7 = clock();    
		// std::cout << "committ update " << ( t7 - t6 ) / (double) CLOCKS_PER_SEC << std::endl;

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
			solve_one_step(i_opt_step, no_timesteps,batch_size,dt,dopt,no_latt_sampling_steps,fname_coll,options);
		};

	};

};