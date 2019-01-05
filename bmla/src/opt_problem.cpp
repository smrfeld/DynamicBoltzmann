#include "../include/bmla_bits/opt_problem.hpp"

// Other headers
#include "../include/bmla_bits/ixn_param.hpp"
#include "../include/bmla_bits/lattice.hpp"
#include "../include/bmla_bits/moment.hpp"
#include "../include/bmla_bits/general.hpp"

#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>

/************************************
* Namespace for bmla
************************************/

namespace bmla {

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

	const std::vector<FName>& FNameColl::get_fnames_all() const {
		return _fnames;
	};
	const FName& FNameColl::get_fname(int idx) const {
		if (idx >= _fnames.size()) {
			std::cerr << ">>> Error: FNameColl::get_fname <<< idx: " << idx << " out of range: " << _fnames.size() << std::endl;
			exit(EXIT_FAILURE);
		};
		return _fnames[idx];
	};

	/********************
	Add fname
	********************/

	void FNameColl::add_fname(FName fname) {
		_fnames.push_back(fname);
		_idxs.push_back(_fnames.size()-1);
	};

	/********************
	Get random subset
	********************/

	std::vector<int> FNameColl::get_random_subset(int size) {
		if (_fnames.size() < size) {
			std::cerr << ">>> Error: FNameColl::get_random_subset <<< size of subset is equal to size of filename collection" << std::endl;
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

	void OptProblem::init_structures(int batch_size) {

		// Ixn params
		for (auto &ixn_param: _ixn_params) {

			// Moment
			ixn_param->get_moment()->set_batch_size(batch_size);
		};
	};

	/********************
	Wake/asleep loop
	********************/

	void OptProblem::wake_sleep_loop(int batch_size, int no_latt_sampling_steps, FNameColl &fname_coll, OptionsWakeSleep options) {
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
				_latt->read_from_file(file.name,file.binary);
			} else {
				// Random lattice...
				_latt->all_units_v_random();
			};

			// Sample hidden
			// Hidden: prob (= false)
			_latt->sample_up_v_to_h(options.layer_wise, options.is_awake_moment_hidden_binary,options.parallel);

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
			for (int i_sampling_step=0; i_sampling_step<no_latt_sampling_steps; i_sampling_step++) 
			{
				// Sample down (hidden -> visible)
				// Visible: binary (= true)
				// Hiddens: binary (= true)
				// _latt->sample_down_h_to_v(layer_wise, true, true);
				_latt->sample_down_h_to_v(options.layer_wise, options.is_asleep_visible_binary, options.is_asleep_hidden_binary,options.parallel);

				// Sample up (visible -> hidden)
				// If not last step:
				// Hiddens: binary (= true)
				// Else:
				// Hiddens: prob (= false)
				if (i_sampling_step != no_latt_sampling_steps-1) {
					_latt->sample_up_v_to_h(options.layer_wise, options.is_asleep_hidden_binary,options.parallel); // binary hiddens
				} else {
					_latt->sample_up_v_to_h(options.layer_wise, options.is_asleep_hidden_binary_final,options.parallel); // prob hiddens
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
	void OptProblem::check_options(int batch_size, double dopt, int no_latt_sampling_steps, OptionsSolve options) {
		// ...
	};

	// One step
	void OptProblem::solve_one_step(int i_opt_step, int batch_size, double dopt, int no_latt_sampling_steps, FNameColl &fname_coll, OptionsSolve options) {

		/*****
		Check options
		*****/

		if (options.should_check_options) {
			check_options(batch_size,dopt,no_latt_sampling_steps,options);
		};

		/*****
		Wake/asleep loop
		*****/

		wake_sleep_loop(batch_size,no_latt_sampling_steps,fname_coll,options.options_wake_sleep);

		if (options.VERBOSE_MOMENT) {
			for (auto &ixn_param: _ixn_params) {
				std::cout << ixn_param->get_name() << " " << std::flush;
				ixn_param->get_moment()->print_moment_comparison();
			};
		};

		/********************
		Form the update
		********************/

		if (options.VERBOSE_UPDATE) {
			std::cout << "--- Calculating update ---" << std::endl;
		};

		double dopt_use = dopt;
		for (auto &ixn_param: _ixn_params) {
			if (!ixn_param->get_is_val_fixed()) {
				// Learning rate
				if (options.MODE_var_learning_rates) {
					dopt_use = options.VAL_var_learning_rates[ixn_param];
				};

				// Update
				if (options.MODE_l2_reg) {
					ixn_param->update_calculate_and_store(dopt_use,options.MODE_l2_reg,options.VAL_l2_lambda[ixn_param],options.VAL_l2_center[ixn_param]);
				} else {
					ixn_param->update_calculate_and_store(dopt_use);
				};
			};
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
			if (!ixn_param->get_is_val_fixed()) {
				ixn_param->update_committ_stored(options.nesterov);
			};
		};

		if (options.VERBOSE_UPDATE) {
			std::cout << "--- [Finished] Committing update ---" << std::endl;
			std::cout << std::endl;
		};
	};

	// Many steps
	void OptProblem::solve(int no_opt_steps, int batch_size, double dopt, int no_latt_sampling_steps, FNameColl &fname_coll, OptionsSolve options) {

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
			solve_one_step(i_opt_step,batch_size,dopt,no_latt_sampling_steps,fname_coll,options);
		};

	};

};