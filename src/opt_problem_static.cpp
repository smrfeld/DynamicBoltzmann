#include "../include/dblz_bits/opt_problem_static.hpp"

// Other headers
#include "../include/dblz_bits/ixn_param.hpp"
#include "../include/dblz_bits/lattice_1d_fully_visible.hpp"
#include "../include/dblz_bits/lattice_centered_hom.hpp"
#include "../include/dblz_bits/moment_diff.hpp"
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
    
    OptProblemStatic::OptProblemStatic() {};
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
    void OptProblemStatic::_move(OptProblemStatic &other) {};
    void OptProblemStatic::_copy(const OptProblemStatic& other) {};
    
    // ***************
    // MARK: - BM PCD
    // ***************
    
    // One step
    void OptProblemStatic::solve_one_step_bm(std::shared_ptr<Lattice> latt, int i_opt_step, int no_steps_awake, int no_steps_asleep, FNameColl &fname_coll, OptionsSolveStatic options, OptionsWakeSleep_BM options_wake_sleep) {
        
        /*****
         Wake/asleep loop
         *****/
        
        int no_markov_awake = latt->get_no_markov_chains(MCType::AWAKE);
        auto fnames = fname_coll.get_random_subset_fnames(no_markov_awake);
        
        latt->wake_sleep_loop_bm(i_opt_step, no_steps_awake, no_steps_asleep, fnames, options_wake_sleep);
        
        // Reap the moments
        latt->reap_ixn_moment_diffs();

        if (options.verbose_moment) {
            for (auto &ixn_param: latt->get_all_ixn_params()) {
                std::cout << ixn_param->get_name() << " " << std::flush;
                ixn_param->get_moment_diff()->print_moment_comparison();
            };
        };
        
        /********************
         Form the update
         ********************/
        
        if (options.verbose_update) {
            std::cout << "--- Calculating update ---" << std::endl;
        };
        
        for (auto &ixn_param: latt->get_all_ixn_params()) {
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
            for (auto &ixn_param: latt->get_all_ixn_params()) {
                if (!ixn_param->get_is_val_fixed()) {
                    ixn_param->update_committ_stored_sgd();
                };
            };
        } else if (options.solver == Solver::ADAM) {
            for (auto &ixn_param: latt->get_all_ixn_params()) {
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
    
    // ***************
    // MARK: - DBM CD centered homogenous params
    // ***************
    
    void OptProblemStatic::solve_one_step_bm(std::shared_ptr<LatticeCenteredHom> lattch, int i_opt_step, int no_steps_awake, int no_steps_asleep, FNameColl &fname_coll, OptionsSolveStatic options, OptionsWakeSleep_BM options_wake_sleep) {
        
        /*****
         Wake/asleep loop
         *****/
        
        int no_markov_awake = lattch->get_no_markov_chains(MCType::AWAKE);
        auto fnames = fname_coll.get_random_subset_fnames(no_markov_awake);
        
        lattch->wake_sleep_loop_bm(i_opt_step, no_steps_awake, no_steps_asleep, fnames, options_wake_sleep);
        
        // Reap the moments
        lattch->reap_ixn_moment_diffs_and_slide_centers(options.sliding_factor);
        
        if (options.verbose_moment) {
            for (auto &ixn_param: lattch->get_all_ixn_params()) {
                std::cout << ixn_param->get_name() << " " << std::flush;
                ixn_param->get_moment_diff()->print_moment_comparison();
            };
        };
        
        /********************
         Form the update
         ********************/
        
        if (options.verbose_update) {
            std::cout << "--- Calculating update ---" << std::endl;
        };
        
        for (auto &ixn_param: lattch->get_all_ixn_params()) {
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
            for (auto &ixn_param: lattch->get_all_ixn_params()) {
                if (!ixn_param->get_is_val_fixed()) {
                    ixn_param->update_committ_stored_sgd();
                };
            };
        } else if (options.solver == Solver::ADAM) {
            for (auto &ixn_param: lattch->get_all_ixn_params()) {
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

    
    // ***************
    // MARK: - RBM CD
    // ***************
    
    void OptProblemStatic::solve_one_step_rbm(std::shared_ptr<Lattice> latt, int i_opt_step, int no_cd_steps, FNameColl &fname_coll, OptionsSolveStatic options, OptionsWakeSleep_RBM options_wake_sleep) {
        
        if (latt->get_no_markov_chains(MCType::AWAKE) != latt->get_no_markov_chains(MCType::ASLEEP)) {
            std::cerr << ">>> OptProblemStatic::solve_one_step_cd <<< Error: no awake and asleep chains must be the same" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        /*****
         Wake/asleep loop
         *****/
        
        int no_markov_awake = latt->get_no_markov_chains(MCType::AWAKE);
        auto fnames = fname_coll.get_random_subset_fnames(no_markov_awake);
        
        latt->wake_sleep_loop_rbm(i_opt_step, no_cd_steps, fnames, options_wake_sleep);
        
        // Reap the moments
        latt->reap_ixn_moment_diffs();

        if (options.verbose_moment) {
            for (auto &ixn_param: latt->get_all_ixn_params()) {
                std::cout << ixn_param->get_name() << " " << std::flush;
                ixn_param->get_moment_diff()->print_moment_comparison();
            };
        };
        
        /********************
         Form the update
         ********************/
        
        if (options.verbose_update) {
            std::cout << "--- Calculating update ---" << std::endl;
        };
        
        for (auto &ixn_param: latt->get_all_ixn_params()) {
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
            for (auto &ixn_param: latt->get_all_ixn_params()) {
                if (!ixn_param->get_is_val_fixed()) {
                    ixn_param->update_committ_stored_sgd();
                };
            };
        } else if (options.solver == Solver::ADAM) {
            for (auto &ixn_param: latt->get_all_ixn_params()) {
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

    // ***************
    // MARK: - RBM CD centered homogenous params
    // ***************
    
    void OptProblemStatic::solve_one_step_rbm(std::shared_ptr<LatticeCenteredHom> lattch, int i_opt_step, int no_cd_steps, FNameColl &fname_coll, OptionsSolveStatic options, OptionsWakeSleep_RBM options_wake_sleep) {
        
        if (lattch->get_no_markov_chains(MCType::AWAKE) != lattch->get_no_markov_chains(MCType::ASLEEP)) {
            std::cerr << ">>> OptProblemStatic::solve_one_step_rbm_centered_hom <<< Error: no awake and asleep chains must be the same" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        /*****
         Wake/asleep loop
         *****/
        
        int no_markov_awake = lattch->get_no_markov_chains(MCType::AWAKE);
        auto fnames = fname_coll.get_random_subset_fnames(no_markov_awake);
        
        lattch->wake_sleep_loop_rbm(i_opt_step, no_cd_steps, fnames, options_wake_sleep);
        
        // Reap the moments
        lattch->reap_ixn_moment_diffs_and_slide_centers(options.sliding_factor);
        
        if (options.verbose_moment) {
            for (auto &ixn_param: lattch->get_all_ixn_params()) {
                std::cout << ixn_param->get_name() << " " << std::flush;
                ixn_param->get_moment_diff()->print_moment_comparison();
            };
        };
        
        /********************
         Form the update
         ********************/
        
        if (options.verbose_update) {
            std::cout << "--- Calculating update ---" << std::endl;
        };
        
        for (auto &ixn_param: lattch->get_all_ixn_params()) {
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
            for (auto &ixn_param: lattch->get_all_ixn_params()) {
                if (!ixn_param->get_is_val_fixed()) {
                    ixn_param->update_committ_stored_sgd();
                };
            };
        } else if (options.solver == Solver::ADAM) {
            for (auto &ixn_param: lattch->get_all_ixn_params()) {
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
    
    // ***************
    // MARK: - 1D Fully visible
    // ***************
    
    void OptProblemStatic::solve_one_step_1d_fully_visible(std::shared_ptr<Lattice1DFullyVisible> latt1dfv, int i_opt_step, int no_cd_steps, FNameColl &fname_coll, OptionsSolveStatic options, OptionsWakeSleep_1DFV_CD options_wake_sleep) {
        
        if (latt1dfv->get_no_markov_chains(MCType::AWAKE) != latt1dfv->get_no_markov_chains(MCType::ASLEEP)) {
            std::cerr << ">>> OptProblemStatic::solve_one_step_1d_fully_visible <<< Error: no awake and asleep chains must be the same" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        /*****
         Wake/asleep loop
         *****/
        
        int no_markov_awake = latt1dfv->get_no_markov_chains(MCType::AWAKE);
        auto fnames = fname_coll.get_random_subset_fnames(no_markov_awake);
        
        latt1dfv->wake_sleep_loop_cd(i_opt_step, no_cd_steps, fnames, options_wake_sleep);
        
        // Reap the moments
        latt1dfv->reap_ixn_moment_diffs();

        if (options.verbose_moment) {
            for (auto &ixn_param: latt1dfv->get_all_ixn_params()) {
                std::cout << ixn_param->get_name() << " " << std::flush;
                ixn_param->get_moment_diff()->print_moment_comparison();
            };
        };
        
        /********************
         Form the update
         ********************/
        
        if (options.verbose_update) {
            std::cout << "--- Calculating update ---" << std::endl;
        };
        
        for (auto &ixn_param: latt1dfv->get_all_ixn_params()) {
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
            for (auto &ixn_param: latt1dfv->get_all_ixn_params()) {
                if (!ixn_param->get_is_val_fixed()) {
                    ixn_param->update_committ_stored_sgd();
                };
            };
        } else if (options.solver == Solver::ADAM) {
            for (auto &ixn_param: latt1dfv->get_all_ixn_params()) {
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
