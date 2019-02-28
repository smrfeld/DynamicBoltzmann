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
    
    OptProblemStatic::OptProblemStatic(std::shared_ptr<Lattice> latt) {
        _latt = latt;
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
        _latt = std::move(other._latt);
    };
    void OptProblemStatic::_copy(const OptProblemStatic& other) {
        if (other._latt) {
            _latt = std::make_shared<Lattice>(*other._latt);
        } else {
            _latt = nullptr;
        };
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
        
        int no_markov_awake = _latt->get_no_markov_chains(MCType::AWAKE);
        auto fnames = fname_coll.get_random_subset_fnames(no_markov_awake);
        
        _latt->wake_sleep_loop(i_opt_step, no_mean_field_updates, no_gibbs_sampling_steps, fnames, options_wake_sleep);
        
        // Reap the moments
        bool slide_means = true;
        if (_latt->get_lattice_mode() == LatticeMode::NORMAL) {
            _latt->reap_moments_and_slide_centers_normal();
        } else if (_latt->get_lattice_mode() == LatticeMode::NORMAL_W_CENTERED_GRADIENT_PT) {
            _latt->reap_moments_and_slide_centers_normal_w_centered_gradient_pt(slide_means,options_wake_sleep.sliding_factor);
        } else if (_latt->get_lattice_mode() == LatticeMode::NORMAL_W_CENTERED_GRADIENT_VEC) {
            _latt->reap_moments_and_slide_centers_normal_w_centered_gradient_vec(slide_means,options_wake_sleep.sliding_factor);
        } else if (_latt->get_lattice_mode() == LatticeMode::CENTERED_PT) {
            _latt->reap_moments_and_slide_centers_centered_pt(slide_means,options_wake_sleep.sliding_factor);
        };
        
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
    
    void OptProblemStatic::solve_one_step_cd(int i_opt_step, int no_cd_steps, FNameColl &fname_coll, OptionsSolveStatic options, OptionsWakeSleep options_wake_sleep) {
        
        if (_latt->get_no_markov_chains(MCType::AWAKE) != _latt->get_no_markov_chains(MCType::ASLEEP)) {
            std::cerr << ">>> OptProblemStatic::solve_one_step_cd <<< Error: no awake and asleep chains must be the same" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        /*****
         Wake/asleep loop
         *****/
        
        int no_markov_awake = _latt->get_no_markov_chains(MCType::AWAKE);
        auto fnames = fname_coll.get_random_subset_fnames(no_markov_awake);
        
        _latt->wake_sleep_loop_cd(i_opt_step, no_cd_steps, fnames, options_wake_sleep);
        
        // Reap the moments
        bool slide_means = true;
        if (_latt->get_lattice_mode() == LatticeMode::NORMAL) {
            _latt->reap_moments_and_slide_centers_normal();
        } else if (_latt->get_lattice_mode() == LatticeMode::NORMAL_W_CENTERED_GRADIENT_PT) {
            _latt->reap_moments_and_slide_centers_normal_w_centered_gradient_pt(slide_means,options_wake_sleep.sliding_factor);
        } else if (_latt->get_lattice_mode() == LatticeMode::NORMAL_W_CENTERED_GRADIENT_VEC) {
            _latt->reap_moments_and_slide_centers_normal_w_centered_gradient_vec(slide_means,options_wake_sleep.sliding_factor);
        } else if (_latt->get_lattice_mode() == LatticeMode::CENTERED_PT) {
            _latt->reap_moments_and_slide_centers_centered_pt(slide_means,options_wake_sleep.sliding_factor);
        };
        
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

};
