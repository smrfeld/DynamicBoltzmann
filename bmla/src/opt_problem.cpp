#include "../include/bmla_bits/opt_problem.hpp"

// Other headers
#include "../include/bmla_bits/ixn_param.hpp"
#include "../include/bmla_bits/lattice.hpp"
#include "../include/bmla_bits/moment.hpp"
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
     OptProblem
     ****************************************/
    
    /********************
     Constructor
     ********************/
    
    OptProblem::OptProblem(std::shared_ptr<Lattice> latt, std::vector<std::shared_ptr<IxnParam>> ixn_params) {
        _latt = latt;
        _ixn_params = ixn_params;
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
        if (other._latt) {
            _latt = std::make_shared<Lattice>(*other._latt);
        };
    };
    
    /********************
     Init structures
     ********************/
    
    void OptProblem::init_structures(int batch_size, int no_markov_chains) {
        
        // Ixn params
        for (auto &ixn_param: _ixn_params) {
            
            // Moment
            ixn_param->get_moment()->set_batch_size(batch_size);
            ixn_param->get_moment()->set_no_markov_chains(no_markov_chains);
        };
        
        // Lattice
        _latt->set_no_markov_chains_asleep(no_markov_chains);
    };
    
    /********************
     Wake/asleep loop
     ********************/
    
    void OptProblem::wake_sleep_loop(int batch_size, int no_markov_chains, int no_cd_sampling_steps, int no_mean_field_updates, FNameColl &fname_coll, OptionsWakeSleep options) {
        if (options.verbose) {
            std::cout << "--- Sampling lattice ---" << std::endl;
        };
        
        // Make a subset
        std::vector<int> idx_subset;
        idx_subset = fname_coll.get_random_subset(batch_size);
        
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
            
            // Switch back to awake stats chain
            _latt->switch_to_markov_chain_awake();
            
            // Read latt
            auto file = fname_coll.get_fname(idx_subset[i_batch]);
            _latt->read_layer_from_file(0,file.name,file.binary);
            
            // Variational inference
            for (auto i=0; i<no_mean_field_updates; i++) {
                _latt->variational_inference_hiddens();
            };
            
            // Reap awake
            _latt->reap_moments(MomentType::AWAKE, i_batch);
            
            // Is START_FROM_DATA, also do the asleep phase
            if (options.cd_mode_asleep == CDModeAsleep::START_FROM_DATA) {
                
                // Start from the current chain
                // Activate up
                _latt->activate_upward_pass(options.options_asleep_start_from_data.start_from_binary_hidden);
                
                // Sample vis, hidden
                for (int i_sampling_step=0; i_sampling_step<no_cd_sampling_steps; i_sampling_step++)
                {
                    if (i_sampling_step != no_cd_sampling_steps-1) {
                        _latt->sample(options.is_asleep_visible_binary, options.is_asleep_hidden_binary);
                    } else {
                        _latt->sample(options.is_asleep_visible_binary, options.is_asleep_hidden_binary_final);
                    };
                };
                
                // Reap asleep
                _latt->reap_moments(MomentType::ASLEEP, i_batch);
                
            };
        };
        
        // If START_FROM_RANDOM or PERSISTENT_CD: run chains
        if (options.cd_mode_asleep == CDModeAsleep::START_FROM_RANDOM) {
            
            // START_FROM_RANDOM
            
            // Run CD sampling
            for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
                
                // Switch to that chain
                _latt->switch_to_markov_chain_asleep(i_chain);
                
                // Random
                _latt->set_random_all_units_in_layer(0,options.options_asleep_start_from_random.start_from_binary_visible);
                _latt->set_random_all_hidden_units(options.options_asleep_start_from_random.start_from_binary_hidden);
                
                // Sample vis, hidden
                for (int i_sampling_step=0; i_sampling_step<no_cd_sampling_steps; i_sampling_step++)
                {
                    if (i_sampling_step != no_cd_sampling_steps-1) {
                        _latt->sample(options.is_asleep_visible_binary, options.is_asleep_hidden_binary);
                    } else {
                        _latt->sample(options.is_asleep_visible_binary, options.is_asleep_hidden_binary_final);
                    };
                };
                
                // Reap asleep
                _latt->reap_moments(MomentType::ASLEEP, i_chain);
            };
            
            // Switch back to awake stats chain
            _latt->switch_to_markov_chain_awake();
            
        } else if (options.cd_mode_asleep == CDModeAsleep::PERSISTENT_CD) {
            
            // PERSISTENT_CD
            
            // Run CD sampling
            for (auto i_chain=0; i_chain<no_markov_chains; i_chain++) {
                
                // Switch to that chain
                _latt->switch_to_markov_chain_asleep(i_chain);
                
                // Sample vis, hidden
                for (int i_sampling_step=0; i_sampling_step<no_cd_sampling_steps; i_sampling_step++)
                {
                    if (i_sampling_step != no_cd_sampling_steps-1) {
                        _latt->sample(options.is_asleep_visible_binary, options.is_asleep_hidden_binary);
                    } else {
                        _latt->sample(options.is_asleep_visible_binary, options.is_asleep_hidden_binary_final);
                    };
                };
                
                // Reap asleep
                _latt->reap_moments(MomentType::ASLEEP, i_chain);
            };
            
            // Switch back to awake stats chain
            _latt->switch_to_markov_chain_awake();
        };
        
        // Average moments
        for (auto &ixn_param: _ixn_params) {
            if (!ixn_param->get_moment()->get_is_awake_moment_fixed()) {
                ixn_param->get_moment()->average_samples(MomentType::AWAKE);
            };
            ixn_param->get_moment()->average_samples(MomentType::ASLEEP);
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
    void OptProblem::check_options(int batch_size, int no_markov_chains, double dopt, int no_cd_sampling_steps, int no_mean_field_updates, OptionsSolve options, OptionsWakeSleep options_wake_sleep) {
        
        if (options_wake_sleep.cd_mode_asleep == CDModeAsleep::START_FROM_DATA) {
            // Must be equal sizes!
            if (batch_size != no_markov_chains) {
                std::cerr << ">>> OptProblemRBM::check_options <<< for CDModeAsleep::START_FROM_DATA, batch size = " << batch_size << " must equal no markov chains = " << no_markov_chains << " because they are initialized at the data." << std::endl;
                exit(EXIT_FAILURE);
            };
        };
    };
    
    // One step
    void OptProblem::solve_one_step(int i_opt_step, int batch_size, int no_markov_chains, double dopt, int no_cd_sampling_steps, int no_mean_field_updates, FNameColl &fname_coll, OptionsSolve options, OptionsWakeSleep options_wake_sleep) {
        
        /*****
         Check options
         *****/
        
        if (options.should_check_options) {
            check_options(batch_size,no_markov_chains,dopt,no_cd_sampling_steps,no_mean_field_updates,options,options_wake_sleep);
        };
        
        /*****
         Wake/asleep loop
         *****/
        
        wake_sleep_loop(batch_size,no_markov_chains,no_cd_sampling_steps,no_mean_field_updates,fname_coll,options_wake_sleep);
        
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
    void OptProblem::solve(int no_opt_steps, int batch_size, int no_markov_chains, double dopt, int no_cd_sampling_steps, int no_mean_field_updates, FNameColl &fname_coll, OptionsSolve options, OptionsWakeSleep options_wake_sleep) {
        
        /*****
         Init structures
         *****/
        
        init_structures(batch_size,no_markov_chains);
        
        /*****
         Go...
         *****/
        
        for (int i_opt_step=1; i_opt_step<=no_opt_steps; i_opt_step++)
        {
            
            std::cout << "------------------" << std::endl;
            std::cout << "Opt step: " << i_opt_step << " / " << no_opt_steps << std::endl;
            std::cout << "------------------" << std::endl;
            
            // Solve
            solve_one_step(i_opt_step,batch_size,no_markov_chains,dopt,no_cd_sampling_steps,no_mean_field_updates,fname_coll,options,options_wake_sleep);
        };
    };
};
