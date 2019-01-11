#include "../include/bmla_bits/opt_problem_rbm.hpp"

// Other headers
#include "../include/bmla_bits/ixn_param.hpp"
#include "../include/bmla_bits/lattice.hpp"
#include "../include/bmla_bits/moment.hpp"
#include "../include/bmla_bits/general.hpp"
#include "../include/bmla_bits/fname.hpp"

#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>

/************************************
* Namespace for bmla
************************************/

namespace bmla {

	/****************************************
	OptProblemRBM
	****************************************/
	
	/********************
	Constructor
	********************/

	OptProblemRBM::OptProblemRBM(std::shared_ptr<Lattice> latt, std::vector<std::shared_ptr<IxnParam>> ixn_params) {
		_latt = latt;
        _ixn_params = ixn_params;
	};
	OptProblemRBM::OptProblemRBM(const OptProblemRBM& other) {
		_copy(other);
	};
	OptProblemRBM::OptProblemRBM(OptProblemRBM&& other) {
		_move(other);
	};
	OptProblemRBM& OptProblemRBM::operator=(const OptProblemRBM& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	OptProblemRBM& OptProblemRBM::operator=(OptProblemRBM&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
	};
	OptProblemRBM::~OptProblemRBM() {
		_clean_up();
	};

	void OptProblemRBM::_clean_up() {};
	void OptProblemRBM::_move(OptProblemRBM &other) {
		_ixn_params = std::move(other._ixn_params);
		_latt = std::move(other._latt);
	};
	void OptProblemRBM::_copy(const OptProblemRBM& other) {
		_ixn_params = other._ixn_params;
		_latt = other._latt;
	};

	/********************
	Init structures
	********************/

	void OptProblemRBM::init_structures(int batch_size) {

		// Ixn params
		for (auto &ixn_param: _ixn_params) {

			// Moment
			ixn_param->get_moment()->set_batch_size(batch_size);
		};
	};

	/********************
	Wake/asleep loop
	********************/

	void OptProblemRBM::wake_sleep_loop(int batch_size, int no_cd_sampling_steps, FNameColl &fname_coll, OptionsWakeSleepRBM options) {
		if (options.verbose) {
			std::cout << "--- Sampling lattice ---" << std::endl;
		};

		// Make a subset
		std::vector<int> idx_subset;
		if (!options.start_with_random_lattice) {
			idx_subset = fname_coll.get_random_subset(batch_size);
		};

		// Reset all moments
		for (auto &ixn_param: _ixn_params) {
			if (!ixn_param->get_moment()->get_is_awake_moment_fixed()) {
				ixn_param->get_moment()->reset_to_zero(MomentType::AWAKE);
			};
			ixn_param->get_moment()->reset_to_zero(MomentType::ASLEEP);
		};

		// Iterate over the batch
		for (int i_batch=0; i_batch<batch_size; i_batch++)
		{

			if (options.verbose) {
				std::cout << "." << std::flush;
			};

			// Read latt
			if (!options.start_with_random_lattice) {
				auto file = fname_coll.get_fname(idx_subset[i_batch]); 
				_latt->read_layer_from_file(0,file.name,file.binary);
			} else {
				// Random lattice...
				_latt->all_units_v_random(true); // binary by default
			};

			// Sample hidden
			// Hidden: prob (= false)
			_latt->sample_rbm_up_v_to_h(options.is_awake_moment_hidden_binary, options.parallel);

			// Reap awake
			for (auto &ixn_param: _ixn_params) {
				if (!ixn_param->get_moment()->get_is_awake_moment_fixed()) {
					ixn_param->get_moment()->reap_in_batch(MomentType::AWAKE, i_batch);
				};
			};

			// Convert hidden units to be binary
			if (!options.is_awake_moment_hidden_binary && options.should_binarize_hidden_after_awake_moment) {
				_latt->all_units_h_binarize();
			};

			// Sample vis, hidden
			for (int i_sampling_step=0; i_sampling_step<no_cd_sampling_steps; i_sampling_step++) 
			{
				// Sample down (hidden -> visible)
				// Visible: binary (= true)
				// Hiddens: binary (= true)
				// _latt->sample_rbm_down_h_to_v(true, true);
				_latt->sample_rbm_down_h_to_v(options.is_asleep_visible_binary,options.parallel);

				// Sample up (visible -> hidden)
				// If not last step:
				// Hiddens: binary (= true)
				// Else:
				// Hiddens: prob (= false)
				if (i_sampling_step != no_cd_sampling_steps-1) {
					_latt->sample_rbm_up_v_to_h(options.is_asleep_hidden_binary,options.parallel); // binary hiddens
				} else {
					_latt->sample_rbm_up_v_to_h(options.is_asleep_hidden_binary_final,options.parallel); // prob hiddens
				};
			};

			// Print
			// latt.print_occupancy();

			// Reap asleep
			for (auto &ixn_param: _ixn_params) {
				// Visible: prob (= false)
				// Hiddens: prob (= false)
				ixn_param->get_moment()->reap_in_batch(MomentType::ASLEEP, i_batch);
			};
		};
		
		// Average moments
		for (auto &ixn_param: _ixn_params) {
			if (!ixn_param->get_moment()->get_is_awake_moment_fixed()) {
				ixn_param->get_moment()->average_reaps(MomentType::AWAKE);
			};
			ixn_param->get_moment()->average_reaps(MomentType::ASLEEP);
		};

		if (options.verbose) {
			std::cout << std::endl;
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
	void OptProblemRBM::check_options(int batch_size, double dopt, int no_cd_sampling_steps, OptionsSolveRBM options) {
		// ...
	};

	// One step
	void OptProblemRBM::solve_one_step(int i_opt_step, int batch_size, double dopt, int no_cd_sampling_steps, FNameColl &fname_coll, OptionsSolveRBM options) {

		/*****
		Check options
		*****/

		if (options.should_check_options) {
			check_options(batch_size,dopt,no_cd_sampling_steps,options);
		};

		/*****
		Wake/asleep loop
		*****/

		wake_sleep_loop(batch_size,no_cd_sampling_steps,fname_coll,options.options_wake_sleep);

		if (options.verbose_moment) {
			for (auto &ixn_param: _ixn_params) {
				std::cout << ixn_param->get_name() << " " << std::flush;
				ixn_param->get_moment()->print_moment_comparison();
			};
		};

		/********************
		Form the update
		********************/

		if (options.verbose_update) {
			std::cout << "--- Calculating update ---" << std::endl;
		};

		for (auto &ixn_param: _ixn_params) {
			if (!ixn_param->get_is_val_fixed()) {
				// Update
				if (options.l2_reg) {
					ixn_param->update_calculate_and_store(options.l2_reg,options.l2_lambda[ixn_param],options.l2_center[ixn_param]);
				} else {
					ixn_param->update_calculate_and_store();
				};
			};
		};

		if (options.verbose_update) {
			std::cout << "--- [Finished] Calculating update ---" << std::endl;
			std::cout << std::endl;
		};

		/********************
		Committ the update
		********************/

		if (options.verbose_update) {
			std::cout << "--- Committing update ---" << std::endl;
		};

		double dopt_use = dopt;
		if (options.solver == Solver::SGD) {
			for (auto &ixn_param: _ixn_params) {
				if (!ixn_param->get_is_val_fixed()) {
					// Learning rate
					if (options.var_learning_rates) {
						dopt_use = options.var_learning_rate_values[ixn_param];
					};

					ixn_param->update_committ_stored_sgd(dopt_use);
				};
			};
		} else if (options.solver == Solver::NESTEROV) {
			for (auto &ixn_param: _ixn_params) {
				if (!ixn_param->get_is_val_fixed()) {
					// Learning rate
					if (options.var_learning_rates) {
						dopt_use = options.var_learning_rate_values[ixn_param];
					};

					ixn_param->update_committ_stored_nesterov(dopt_use,options.nesterov_acc);
				};
			};
		} else if (options.solver == Solver::ADAM) {
			for (auto &ixn_param: _ixn_params) {
				if (!ixn_param->get_is_val_fixed()) {
					// Learning rate
					if (options.var_learning_rates) {
						dopt_use = options.var_learning_rate_values[ixn_param];
					};

					ixn_param->update_committ_stored_adam(dopt_use,i_opt_step,options.adam_beta_1,options.adam_beta_2,options.adam_eps);
				};
			};			
		};

		if (options.verbose_update) {
			std::cout << "--- [Finished] Committing update ---" << std::endl;
			std::cout << std::endl;
		};
	};

	// Many steps
	void OptProblemRBM::solve(int no_opt_steps, int batch_size, double dopt, int no_cd_sampling_steps, FNameColl &fname_coll, OptionsSolveRBM options) {

		/*****
		Init structures
		*****/

		init_structures(batch_size);

		/*****
		Go...
		*****/

		for (int i_opt_step=1; i_opt_step<=no_opt_steps; i_opt_step++)
		{

			std::cout << "------------------" << std::endl;
			std::cout << "Opt step: " << i_opt_step << " / " << no_opt_steps << std::endl;
			std::cout << "------------------" << std::endl;

			// Solve
			solve_one_step(i_opt_step,batch_size,dopt,no_cd_sampling_steps,fname_coll,options);
		};

	};

};
