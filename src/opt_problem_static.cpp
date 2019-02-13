#include "../include/dblz_bits/opt_problem_static.hpp"

// Other headers
#include "../include/dblz_bits/ixn_param.hpp"
#include "../include/dblz_bits/lattice.hpp"
#include "../include/dblz_bits/moment.hpp"
#include "../include/dblz_bits/fname.hpp"
#include "../include/dblz_bits/general.hpp"

#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>

/************************************
 * Namespace for bmla
 ************************************/

namespace dblz {
    
    /****************************************
     OptProblemStatic
     ****************************************/
    
    /********************
     Constructor
     ********************/
    
    OptProblemStatic::OptProblemStatic(std::shared_ptr<Lattice> latt, int no_markov_chains_awake, int no_markov_chains_asleep) {
        _latt = latt;
        
        set_no_markov_chains(MCType::AWAKE,no_markov_chains_awake);
        set_no_markov_chains(MCType::ASLEEP,no_markov_chains_asleep);
    };
    OptProblemStatic::OptProblemStatic(const OptProblemStatic& other) {
        _copy(other);
    };
    OptProblemStatic::OptProblemStatic(OptProblemStatic&& other) {
        _move(other);
    };
    OptProblemStatic& OptProblemStatic::operator=(const OptProblemStatic& other) {
        if (this != &other) {
            _clean_up();
            _copy(other);
        };
        return *this;
    };
    OptProblemStatic& OptProblemStatic::operator=(OptProblemStatic&& other) {
        if (this != &other) {
            _clean_up();
            _move(other);
        };
        return *this;
    };
    OptProblemStatic::~OptProblemStatic() {
        _clean_up();
    };
    
    void OptProblemStatic::_clean_up() {};
    void OptProblemStatic::_move(OptProblemStatic &other) {
        _no_markov_chains = other._no_markov_chains;
        _latt = std::move(other._latt);
    };
    void OptProblemStatic::_copy(const OptProblemStatic& other) {
        _no_markov_chains = other._no_markov_chains;
        if (other._latt) {
            _latt = std::make_shared<Lattice>(*other._latt);
        };
    };
    
    /********************
     Init structures
     ********************/

    void OptProblemStatic::set_no_markov_chains(MCType chain, int no_markov_chains) {
        
        // Store
        _no_markov_chains[chain] = no_markov_chains;
        
        // Lattice
        _latt->set_no_markov_chains(chain, no_markov_chains);
    };
    
    /********************
     Solve
     ********************/
    
    // Check if options passed are valid
    void OptProblemStatic::check_options(int no_mean_field_updates, int no_gibbs_sampling_steps, OptionsSolveStatic options, OptionsWakeSleep options_wake_sleep) {
    };
    
    // One step
    void OptProblemStatic::solve_one_step(int i_opt_step, int no_mean_field_updates, int no_gibbs_sampling_steps, FNameColl &fname_coll, OptionsSolveStatic options, OptionsWakeSleep options_wake_sleep) {
        
        /*****
         Check options
         *****/
        
        if (options.should_check_options) {
            check_options(no_mean_field_updates,no_gibbs_sampling_steps,options,options_wake_sleep);
        };
        
        /*****
         Wake/asleep loop
         *****/
        
        auto fnames = fname_coll.get_random_subset_fnames(_no_markov_chains.at(MCType::AWAKE));
        
        _latt->wake_sleep_loop(i_opt_step, no_mean_field_updates, no_gibbs_sampling_steps, fnames, options_wake_sleep);
        
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
                    ixn_param->update_calculate_and_store_l2(options.l2_lambda[ixn_param],options.l2_center[ixn_param]);
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
        
        if (options.solver == Solver::SGD) {
            for (auto &ixn_param: _latt->get_all_ixn_params()) {
                if (!ixn_param->get_is_val_fixed()) {
                    ixn_param->update_committ_stored_sgd();
                };
            };
        } else if (options.solver == Solver::NESTEROV) {
            for (auto &ixn_param: _latt->get_all_ixn_params()) {
                if (!ixn_param->get_is_val_fixed()) {
                    ixn_param->update_committ_stored_nesterov(options.nesterov_acc);
                };
            };
        } else if (options.solver == Solver::ADAM) {
            for (auto &ixn_param: _latt->get_all_ixn_params()) {
                if (!ixn_param->get_is_val_fixed()) {
                    ixn_param->update_committ_stored_adam(i_opt_step,options.adam_beta_1,options.adam_beta_2,options.adam_eps);
                };
            };
        };
        
        if (options.verbose_update) {
            std::cout << "--- [Finished] Committing update ---" << std::endl;
            std::cout << std::endl;
        };
    };
    
    // Many steps
    void OptProblemStatic::solve(int no_opt_steps, int no_mean_field_updates, int no_gibbs_sampling_steps, FNameColl &fname_coll, OptionsSolveStatic options, OptionsWakeSleep options_wake_sleep) {
        
        for (int i_opt_step=1; i_opt_step<=no_opt_steps; i_opt_step++)
        {
            
            std::cout << "------------------" << std::endl;
            std::cout << "Opt step: " << i_opt_step << " / " << no_opt_steps << std::endl;
            std::cout << "------------------" << std::endl;
            
            // Solve
            solve_one_step(i_opt_step,no_mean_field_updates,no_gibbs_sampling_steps,fname_coll,options,options_wake_sleep);
        };
    };
};
