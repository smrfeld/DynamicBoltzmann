#include "../include/bmla_bits/opt_problem.hpp"

// Other headers
#include "../include/bmla_bits/ixn_param.hpp"
#include "../include/bmla_bits/lattice.hpp"
#include "../include/bmla_bits/moment.hpp"
#include "../include/bmla_bits/fname.hpp"
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
     OptProblem
     ****************************************/
    
    /********************
     Constructor
     ********************/
    
    OptProblem::OptProblem(std::shared_ptr<Lattice> latt, int no_markov_chains_awake, int no_markov_chains_asleep) {
        _latt = latt;
        
        set_no_markov_chains(MCType::AWAKE,no_markov_chains_awake);
        set_no_markov_chains(MCType::ASLEEP,no_markov_chains_asleep);
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
        _no_markov_chains = other._no_markov_chains;
        _latt = std::move(other._latt);
    };
    void OptProblem::_copy(const OptProblem& other) {
        _no_markov_chains = other._no_markov_chains;
        if (other._latt) {
            _latt = std::make_shared<Lattice>(*other._latt);
        };
    };
    
    /********************
     Init structures
     ********************/

    void OptProblem::set_no_markov_chains(MCType chain, int no_markov_chains) {
        
        // Store
        _no_markov_chains[chain] = no_markov_chains;
        
        // Moments
        for (auto &ixn_param: _latt->get_all_ixn_params()) {
            ixn_param->get_moment()->set_no_markov_chains(chain, no_markov_chains);
        };
        
        // Lattice
        _latt->set_no_markov_chains(chain, no_markov_chains);
    };
    
    /********************
     Wake/asleep loop
     ********************/
    
    void OptProblem::wake_sleep_loop(int i_opt_step, int no_mean_field_updates, int no_gibbs_sampling_steps, FNameColl &fname_coll, OptionsWakeSleep options) {
        if (options.verbose) {
            std::cout << "--- Sampling lattice ---" << std::endl;
        };
        
        // Reset all moments
        for (auto &ixn_param: _latt->get_all_ixn_params()) {
            if (!ixn_param->get_moment()->get_is_awake_moment_fixed()) {
                ixn_param->get_moment()->reset_to_zero(MCType::AWAKE);
            };
            ixn_param->get_moment()->reset_to_zero(MCType::ASLEEP);
        };
        
        // AWAKE PHASE
        
        // Make a batch subset
        std::vector<int> idx_subset;
        idx_subset = fname_coll.get_random_subset(_no_markov_chains[MCType::AWAKE]);
        
        // Read in the batch
        for (int i_chain=0; i_chain<_no_markov_chains[MCType::AWAKE]; i_chain++)
        {
            auto file = fname_coll.get_fname(idx_subset[i_chain]);
            _latt->read_layer_from_file(MCType::AWAKE, i_chain, 0, file.name, file.binary);
            
        };
        
        // Option (1): init MF with random hidden layers with prob units
        /*
        for (int i_chain=0; i_chain<_no_markov_chains[MCType::AWAKE]; i_chain++) {
            _latt->set_random_all_hidden_units(MCType::AWAKE, i_chain, false);
        };
         */
        // Option (2): upward pass with 2x weights (DBM) to activate probabilitsic units
        // (faster to converge!!!)
        _latt->activate_upward_pass_with_2x_weights_1x_bias(MCType::AWAKE, false);
        
        // Variational inference
        for (auto i=0; i<no_mean_field_updates; i++) {
            _latt->mean_field_hiddens_step();
        };
            
        // Reap awake
        _latt->reap_moments(MCType::AWAKE);
        
        // Write out the lattices
        if (options.write_after_awake) {
            for (auto i_chain=0; i_chain<_no_markov_chains[MCType::AWAKE]; i_chain++) {
                for (auto layer=0; layer<_latt->get_no_layers(); layer++) {
                    bool write_binary = false;
                    if (layer == 0) {
                        write_binary = true;
                    };
                    _latt->write_layer_to_file(MCType::AWAKE, i_chain, layer, options.write_after_awake_dir+"/"+pad_str(i_chain,2)+"_"+pad_str(layer,2)+".txt", write_binary);
                };
            };
        };
        
        // ASLEEP PHASE - PERSISTENT_CD
        
        // Run CD sampling
        for (auto i_chain=0; i_chain<_no_markov_chains[MCType::ASLEEP]; i_chain++) {
            
            // Sample vis, hidden
            for (int i_sampling_step=0; i_sampling_step<no_gibbs_sampling_steps; i_sampling_step++)
            {
                if (i_sampling_step != no_gibbs_sampling_steps-1) {
                    _latt->gibbs_sampling_step(options.is_asleep_visible_binary, options.is_asleep_hidden_binary);
                } else {
                    if (options.is_asleep_visible_binary_final && options.is_asleep_hidden_binary_final) {
                        // All binary
                        _latt->gibbs_sampling_step(options.is_asleep_visible_binary_final, options.is_asleep_hidden_binary_final);
                    } else {
                        // Parallel for non-binary options
                        _latt->gibbs_sampling_step_parallel(options.is_asleep_visible_binary_final, options.is_asleep_hidden_binary_final);
                    };
                };
            };
        };
        
        // Reap asleep
        _latt->reap_moments(MCType::ASLEEP);

        // Write out the lattices
        if (options.write_after_asleep) {
            for (auto i_chain=0; i_chain<_no_markov_chains[MCType::ASLEEP]; i_chain++) {
                for (auto layer=0; layer<_latt->get_no_layers(); layer++) {
                    bool write_binary = false;
                    if (layer == 0) {
                        write_binary = true;
                    };
                    _latt->write_layer_to_file(MCType::ASLEEP, i_chain, layer, options.write_after_asleep_dir+"/"+pad_str(i_chain,2)+"_"+pad_str(layer,2)+".txt", write_binary);
                };
            };
        };
        
        // Average moments
        for (auto &ixn_param: _latt->get_all_ixn_params()) {
            // std::cout << "averaging: " << ixn_param->get_name() << std::endl;
            if (!ixn_param->get_moment()->get_is_awake_moment_fixed()) {
                ixn_param->get_moment()->average_moment_samples(MCType::AWAKE);
            };
            ixn_param->get_moment()->average_moment_samples(MCType::ASLEEP);
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
    void OptProblem::check_options(double dopt, int no_mean_field_updates, int no_gibbs_sampling_steps, OptionsSolve options, OptionsWakeSleep options_wake_sleep) {
    };
    
    // One step
    void OptProblem::solve_one_step(int i_opt_step, double dopt, int no_mean_field_updates, int no_gibbs_sampling_steps, FNameColl &fname_coll, OptionsSolve options, OptionsWakeSleep options_wake_sleep) {
        
        /*****
         Check options
         *****/
        
        if (options.should_check_options) {
            check_options(dopt,no_mean_field_updates,no_gibbs_sampling_steps,options,options_wake_sleep);
        };
        
        /*****
         Wake/asleep loop
         *****/
        
        wake_sleep_loop(i_opt_step,no_mean_field_updates,no_gibbs_sampling_steps,fname_coll,options_wake_sleep);
        
        if (options.verbose_moment) {
            for (auto &ixn_param: _latt->get_all_ixn_params()) {
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
        
        for (auto &ixn_param: _latt->get_all_ixn_params()) {
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
            for (auto &ixn_param: _latt->get_all_ixn_params()) {
                if (!ixn_param->get_is_val_fixed()) {
                    // Learning rate
                    if (options.var_learning_rates) {
                        dopt_use = options.var_learning_rate_values[ixn_param];
                    };
                    
                    ixn_param->update_committ_stored_sgd(dopt_use);
                };
            };
        } else if (options.solver == Solver::NESTEROV) {
            for (auto &ixn_param: _latt->get_all_ixn_params()) {
                if (!ixn_param->get_is_val_fixed()) {
                    // Learning rate
                    if (options.var_learning_rates) {
                        dopt_use = options.var_learning_rate_values[ixn_param];
                    };
                    
                    ixn_param->update_committ_stored_nesterov(dopt_use,options.nesterov_acc);
                };
            };
        } else if (options.solver == Solver::ADAM) {
            for (auto &ixn_param: _latt->get_all_ixn_params()) {
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
    void OptProblem::solve(int no_opt_steps, double dopt, int no_mean_field_updates, int no_gibbs_sampling_steps, FNameColl &fname_coll, OptionsSolve options, OptionsWakeSleep options_wake_sleep) {
        
        for (int i_opt_step=1; i_opt_step<=no_opt_steps; i_opt_step++)
        {
            
            std::cout << "------------------" << std::endl;
            std::cout << "Opt step: " << i_opt_step << " / " << no_opt_steps << std::endl;
            std::cout << "------------------" << std::endl;
            
            // Solve
            solve_one_step(i_opt_step,dopt,no_mean_field_updates,no_gibbs_sampling_steps,fname_coll,options,options_wake_sleep);
        };
    };
};
