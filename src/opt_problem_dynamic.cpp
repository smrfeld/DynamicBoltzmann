#include "../include/dblz_bits/opt_problem_dynamic.hpp"

// Other headers
#include "../include/dblz_bits/ixn_param.hpp"
#include "../include/dblz_bits/ixn_param_traj.hpp"
#include "../include/dblz_bits/lattice_1d_fully_visible.hpp"
#include "../include/dblz_bits/lattice_traj_1d_fully_visible.hpp"
#include "../include/dblz_bits/lattice_alternating_binary.hpp"
#include "../include/dblz_bits/lattice_traj_alternating_binary.hpp"
#include "../include/dblz_bits/lattice_centered_hom.hpp"
#include "../include/dblz_bits/lattice_traj_centered_hom.hpp"
#include "../include/dblz_bits/moment_diff.hpp"
#include "../include/dblz_bits/fname_traj.hpp"
#include "../include/dblz_bits/general.hpp"
#include "../include/dblz_bits/adjoint_obs.hpp"
#include "../include/dblz_bits/adjoint_params.hpp"
#include "../include/dblz_bits/adjoint_params_centered_hom.hpp"
#include "../include/dblz_bits/diff_eq_rhs.hpp"
#include "../include/dblz_bits/species.hpp"

#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <ctime>

/************************************
 * Namespace for bmla
 ************************************/

namespace dblz {
    
    /****************************************
     OptProblemDynamic
     ****************************************/
    
    /********************
     Constructor
     ********************/
    
    OptProblemDynamic::OptProblemDynamic() {};
    OptProblemDynamic::OptProblemDynamic(const OptProblemDynamic& other) {
        _copy(other);
    };
    OptProblemDynamic::OptProblemDynamic(OptProblemDynamic&& other) {
        _move(other);
    };
    OptProblemDynamic& OptProblemDynamic::operator=(const OptProblemDynamic& other) {
        if (this != &other) {
            _clean_up();
            _copy(other);
        };
        return *this;
    };
    OptProblemDynamic& OptProblemDynamic::operator=(OptProblemDynamic&& other) {
        if (this != &other) {
            _clean_up();
            _move(other);
        };
        return *this;
    };
    OptProblemDynamic::~OptProblemDynamic() {
        _clean_up();
    };
    
    void OptProblemDynamic::_clean_up() {};
    void OptProblemDynamic::_move(OptProblemDynamic &other) {};
    void OptProblemDynamic::_copy(const OptProblemDynamic& other) {};
    
    // ***************
    // MARK: - Solve ixn params
    // ***************
    
    void OptProblemDynamic::solve_ixn_param_trajs_step(const std::vector<std::shared_ptr<IxnParamTraj>> &ixn_param_trajs, const std::vector<std::shared_ptr<Domain>> &domains, double dt, int timepoint) const {
        
        // clock_t t0 = clock();

        // Form abscissas
        for (auto domain: domains) {
            domain->calculate_val_at_timepoint(timepoint);
        };
        
        // clock_t t1 = clock();

        // Solve
        for (auto ixn_param_traj: ixn_param_trajs) {
            if (!ixn_param_traj->get_is_val_fixed()) {
                ixn_param_traj->solve_diff_eq_at_timepoint_to_plus_one(timepoint, dt,false);
            };
        };
        
        // clock_t t2 = clock();
        /*
        double dt1 = (t1-t0)  / (double) CLOCKS_PER_SEC;
        double dt2 = (t2-t1)  / (double) CLOCKS_PER_SEC;
        double dt_tot = dt1 + dt2;
        std::cout << "[time " << dt_tot << "] [ " << dt1/dt_tot << "] [ " << dt2/dt_tot << "]" << std::endl;
         */
    };
    void OptProblemDynamic::solve_ixn_param_trajs_step(const std::vector<std::shared_ptr<IxnParamTraj>> &ixn_param_trajs, const std::vector<std::shared_ptr<Domain>> &domains, double dt, int timepoint, int no_steps_per_step) const {
        
        // clock_t t0 = clock();
        
        // Get initial vals at this timepoint
        for (auto ixn_param_traj: ixn_param_trajs) {
            ixn_param_traj->set_substep_val_current(ixn_param_traj->get_ixn_param_at_timepoint(timepoint)->get_val());
        };
        
        // clock_t t1 = clock();
        
        // Calculate diff eqs to a new step
        for (auto i=0; i<no_steps_per_step; i++) {
            
            // Calculate domain vals
            for (auto domain: domains) {
                domain->calculate_substep_val();
            };
            
            // New substep vals
            for (auto ixn_param_traj: ixn_param_trajs) {
                if (!ixn_param_traj->get_is_val_fixed()) {
                    // std::cout << ixn_param_traj->get_name() << " " << timepoint << " " << i << " " << ixn_param_traj->get_diff_eq_rhs()->get_val_from_map(vals, i) << std::endl;
                    ixn_param_traj->set_substep_val_new(ixn_param_traj->get_substep_val() + (dt / no_steps_per_step) * ixn_param_traj->get_diff_eq_rhs()->get_substep_val(false));
                };
            };
            
            // Committ
            for (auto ixn_param_traj: ixn_param_trajs) {
                if (!ixn_param_traj->get_is_val_fixed()) {
                    ixn_param_traj->move_to_new_substep_val();
                };
            };
        };
        
        // clock_t t2 = clock();
        
        // Write final vals
        for (auto ixn_param_traj: ixn_param_trajs) {
            if (!ixn_param_traj->get_is_val_fixed()) {                ixn_param_traj->get_ixn_param_at_timepoint(timepoint+1)->set_val(ixn_param_traj->get_substep_val());
            };
        };
        
        // clock_t t3 = clock();
        
        /*
        double dt1 = (t1-t0)  / (double) CLOCKS_PER_SEC;
        double dt2 = (t2-t1)  / (double) CLOCKS_PER_SEC;
        double dt3 = (t3-t2)  / (double) CLOCKS_PER_SEC;
        double dt_tot = dt1 + dt2 + dt3;
        std::cout << "[time " << dt_tot << "] [ " << dt1/dt_tot << "] [ " << dt2/dt_tot << "] [ " << dt3/dt_tot << "]" << std::endl;
         */
    };

    void OptProblemDynamic::solve_ixn_param_trajs(const std::vector<std::shared_ptr<IxnParamTraj>> &ixn_param_trajs, const std::vector<std::shared_ptr<Domain>> &domains, double dt, int timepoint_start, int no_timesteps) const {
        
        for (auto timepoint=timepoint_start; timepoint<timepoint_start+no_timesteps; timepoint++) {
            solve_ixn_param_trajs_step(ixn_param_trajs, domains, dt, timepoint);
        };
        
        // Set domain val at final timestep
        for (auto domain: domains) {
            domain->calculate_val_at_timepoint(timepoint_start+no_timesteps);
        };
        
        // Now all the domain vals have been calculated!
    };
    void OptProblemDynamic::solve_ixn_param_trajs(const std::vector<std::shared_ptr<IxnParamTraj>> &ixn_param_trajs, const std::vector<std::shared_ptr<Domain>> &domains, double dt, int timepoint_start, int no_timesteps, int no_steps_per_step) const {
        
        for (auto timepoint=timepoint_start; timepoint<timepoint_start+no_timesteps; timepoint++) {
            // Set domain val
            for (auto domain: domains) {
                domain->calculate_val_at_timepoint(timepoint);
            };

            // Solve forward
            solve_ixn_param_trajs_step(ixn_param_trajs, domains, dt, timepoint, no_steps_per_step);
        };
        
        // Set domain val at final timestep
        for (auto domain: domains) {
            domain->calculate_val_at_timepoint(timepoint_start+no_timesteps);
        };
        
        // Now all the domain vals have been calculated!
    };

    // ***************
    // MARK: - BM PCD params
    // ***************

    void OptProblemDynamic::solve_one_step_bm_params(std::shared_ptr<LatticeTrajCenteredHom> latt_traj, const std::vector<std::shared_ptr<IxnParamTraj>> &all_ixn_param_trajs, const std::vector<std::shared_ptr<Domain>> &all_domains, int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_steps_awake, int no_steps_asleep, FNameTrajColl &fname_traj_coll, const std::vector<std::shared_ptr<AdjointParamsCenteredHomDerivTerm>> &all_deriv_terms, OptionsSolveDynamic options, OptionsWakeSleep_BM options_wake_sleep) {

    solve_one_step_bm_params_without_committ(latt_traj,all_ixn_param_trajs,all_domains,i_opt_step,timepoint_start_SIP,no_timesteps_SIP,timepoint_start_WS,no_timesteps_WS,timepoint_start_A,no_timesteps_A,dt,no_steps_awake,no_steps_asleep,fname_traj_coll,all_deriv_terms,options,options_wake_sleep);
        
        committ_step(latt_traj->get_all_ixn_param_trajs(), i_opt_step, options);
    };
    
    
    void OptProblemDynamic::solve_one_step_bm_params_without_committ(std::shared_ptr<LatticeTrajCenteredHom> latt_traj, const std::vector<std::shared_ptr<IxnParamTraj>> &all_ixn_param_trajs, const std::vector<std::shared_ptr<Domain>> &all_domains, int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_steps_awake, int no_steps_asleep, FNameTrajColl &fname_traj_coll, const std::vector<std::shared_ptr<AdjointParamsCenteredHomDerivTerm>> &all_deriv_terms, OptionsSolveDynamic options, OptionsWakeSleep_BM options_wake_sleep) {
        
        if (options.locking_mode) {
            std::cerr << ">>> OptProblemDynamic::solve_one_step_bm_params_without_committ <<< Locking mode not supported here" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        if (options.l2_reg || options.l2_reg_traj) {
            std::cerr << ">>> OptProblemDynamic::solve_one_step_bm_params_without_committ <<< L2 reg mode not supported here" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        /*****
         Solve diff eq for F
         *****/
        
        clock_t t0 = clock();
        
        if (options.no_steps_per_step_IP == 1) {
            solve_ixn_param_trajs(all_ixn_param_trajs, all_domains, dt, timepoint_start_SIP, no_timesteps_SIP);
        } else {
            solve_ixn_param_trajs(all_ixn_param_trajs, all_domains, dt, timepoint_start_SIP, no_timesteps_SIP, options.no_steps_per_step_IP);
        };

        // After now, no more need to reform coordinates
        bool form_abscissas = false;
        
        /*****
         Wake/asleep loop
         *****/
        
        clock_t t1 = clock();
        
        int no_awake_chains = latt_traj->get_no_markov_chains(MCType::AWAKE);
        std::vector<std::vector<FName>> fname_coll = fname_traj_coll.get_random_subset_fnames(no_awake_chains, timepoint_start_WS, no_timesteps_WS);

        for (auto timepoint=timepoint_start_WS; timepoint<=timepoint_start_WS+no_timesteps_WS; timepoint++) {
            latt_traj->get_lattice_at_timepoint(timepoint)->wake_sleep_loop_bm(i_opt_step, no_steps_awake, no_steps_asleep, fname_coll.at(timepoint-timepoint_start_WS), options_wake_sleep);
        };
        
        clock_t t2 = clock();
        
        // Reap the moments
        bool calculate_offset = false;
        for (auto timepoint=timepoint_start_WS; timepoint<=timepoint_start_WS+no_timesteps_WS; timepoint++) {
            latt_traj->get_lattice_centered_hom_at_timepoint(timepoint)->reap_ixn_moment_diffs_and_slide_centers(options.sliding_factor,calculate_offset);
        };
        
        // Print
        if (options.verbose) {
            for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
                std::cout << ixn_param_traj->get_name() << " moments [" << timepoint_start_WS << "," << timepoint_start_WS+no_timesteps_WS << "]" << std::endl;
                ixn_param_traj->print_moment_diff_traj(timepoint_start_WS, no_timesteps_WS);
            };
        };

        /********************
         Solve diff eq for adjoint
         ********************/
        
        clock_t t3 = clock();
        
        // Set zero endpoint
        for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
            ixn_param_traj->get_adjoint()->set_timepoint_zero_end_cond(timepoint_start_A + no_timesteps_A);
        };

        for (auto timepoint=timepoint_start_A + no_timesteps_A; timepoint>timepoint_start_A; timepoint--) {
            // Calculate deriv terms
            for (auto deriv_term: all_deriv_terms) {
                deriv_term->calculate_val_at_timepoint(timepoint,form_abscissas);
            };
                        
            // Solve diff eq
            for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
                if (!ixn_param_traj->get_is_val_fixed()) {
                    ixn_param_traj->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one(timepoint, dt, form_abscissas);
                };
            };
        };
        
        /********************
         Form the update
         ********************/
        
        clock_t t4 = clock();
        
        for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
            if (!ixn_param_traj->get_is_val_fixed()) {
                ixn_param_traj->get_diff_eq_rhs()->update_calculate_and_store(timepoint_start_A,timepoint_start_A+no_timesteps_A,dt, form_abscissas);
            };
        };
        
        clock_t t5 = clock();

        if (options.verbose_timing) {
            double dt1 = (t1-t0)  / (double) CLOCKS_PER_SEC;
            double dt2 = (t2-t1)  / (double) CLOCKS_PER_SEC;
            double dt3 = (t3-t2)  / (double) CLOCKS_PER_SEC;
            double dt4 = (t4-t3)  / (double) CLOCKS_PER_SEC;
            double dt5 = (t5-t4)  / (double) CLOCKS_PER_SEC;
            double dt_tot = dt1 + dt2 + dt3 + dt4 + dt5;
            std::cout << "[time " << dt_tot << "] [F " << dt1/dt_tot << "] [wake/sleep " << dt2/dt_tot << "] [reap " << dt3/dt_tot << "] [adj " << dt4/dt_tot << "] [form update " << dt5/dt_tot << "]" << std::endl;
        };
    };
    
    // ***************
    // MARK: - RBM CD params
    // ***************
    
    void OptProblemDynamic::solve_one_step_rbm_params(std::shared_ptr<LatticeTraj> latt_traj, const std::vector<std::shared_ptr<IxnParamTraj>> &all_ixn_param_trajs, const std::vector<std::shared_ptr<Domain>> &all_domains, int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_cd_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_RBM options_wake_sleep) {
        
        solve_one_step_rbm_params_without_committ(latt_traj,all_ixn_param_trajs,all_domains,i_opt_step,timepoint_start_SIP,no_timesteps_SIP,timepoint_start_WS,no_timesteps_WS,timepoint_start_A,no_timesteps_A,dt,no_cd_steps,fname_traj_coll,options,options_wake_sleep);
        
        committ_step(latt_traj->get_all_ixn_param_trajs(), i_opt_step, options);
    };
    
    
    void OptProblemDynamic::solve_one_step_rbm_params_without_committ(std::shared_ptr<LatticeTraj> latt_traj, const std::vector<std::shared_ptr<IxnParamTraj>> &all_ixn_param_trajs, const std::vector<std::shared_ptr<Domain>> &all_domains, int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_cd_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_RBM options_wake_sleep) {
        
        if (latt_traj->get_no_markov_chains(MCType::AWAKE) != latt_traj->get_no_markov_chains(MCType::ASLEEP)) {
            std::cerr << ">>> OptProblemDynamic::solve_one_step_rbm_params_without_committ <<< asleep chains must equal awake chains" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        /*****
         Solve diff eq for F
         *****/
        
        clock_t t0 = clock();
        
        solve_ixn_param_trajs(all_ixn_param_trajs,all_domains, dt, timepoint_start_SIP, no_timesteps_SIP, options.no_steps_per_step_IP);

        // After now, no more need to reform coordinates
        bool form_abscissas = false;
        
        /*****
         Adjust adjoint range
         *****/
        
        int timepoint_start_A_use = timepoint_start_A;
        int no_timesteps_A_use = no_timesteps_A;
        int timepoint_start_WS_use = timepoint_start_WS;
        int no_timesteps_WS_use = no_timesteps_WS;
        
        if (options.locking_mode) {
            
            // Look at the cell in timepoint_start_A_use
            // Extend adjoint range to any timesteps also in this cell
            q3c1::Cell* cell = nullptr, *cell_prev = nullptr;
            bool repeat = false;
            if (timepoint_start_A_use != 0) {
                do {
                    // Dont repeat
                    repeat = false;
                    // Run through all ixn params
                    for (auto ixn: latt_traj->get_all_ixn_param_trajs()) {
                        // Get the cell at timepoint_start_A_use
                        cell = ixn->get_diff_eq_rhs()->get_cell_at_timepoint(timepoint_start_A_use);
                        cell_prev = ixn->get_diff_eq_rhs()->get_cell_at_timepoint(timepoint_start_A_use-1);
                        if (cell == cell_prev) {
                            // Same cell, extend adjoint domain
                            timepoint_start_A_use -= 1;
                            no_timesteps_A_use += 1;
                            // Ensure that this is covered by the wake-sleep range
                            if (timepoint_start_A_use < timepoint_start_WS_use) {
                                timepoint_start_WS_use = timepoint_start_A_use;
                            };
                            if (timepoint_start_WS_use + no_timesteps_WS_use < timepoint_start_A_use + no_timesteps_A_use) {
                                no_timesteps_WS_use = timepoint_start_A_use + no_timesteps_A_use - timepoint_start_WS_use;
                            };
                            // Restart the search
                            repeat = true;
                            break;
                        };
                    };
                } while (repeat && timepoint_start_A_use != 0);
            };
            std::cout << "Adjoint range is: " << timepoint_start_A_use << " to: " << timepoint_start_A_use + no_timesteps_A_use << std::endl;
            std::cout << "Wake/sleep range is: " << timepoint_start_WS_use << " to: " << timepoint_start_WS_use + no_timesteps_WS_use << std::endl;
            
            // Free all cells at timepoint start and after
            // Note: order matters, do this first!
            for (auto timepoint=timepoint_start_A_use; timepoint<=timepoint_start_A_use+no_timesteps_A_use; timepoint++) {
                for (auto ixn: latt_traj->get_all_ixn_param_trajs()) {
                    ixn->get_diff_eq_rhs()->fix_all_verts_around_at_timepoint(timepoint, false);
                };
            };

            // Fix all cells before the adjoint timepoint start
            // Note: order matters, do this second!
            for (auto timepoint=0; timepoint<timepoint_start_A_use; timepoint++) {
                for (auto ixn: latt_traj->get_all_ixn_param_trajs()) {
                    ixn->get_diff_eq_rhs()->fix_all_verts_around_at_timepoint(timepoint, true);
                };
            };
        };
        
        /*****
         Wake/asleep loop
         *****/
        
        clock_t t1 = clock();
        
        int no_awake_chains = latt_traj->get_no_markov_chains(MCType::AWAKE);
        std::vector<std::vector<FName>> fname_coll = fname_traj_coll.get_random_subset_fnames(no_awake_chains, timepoint_start_WS_use, no_timesteps_WS_use);
        
        for (auto timepoint=timepoint_start_WS_use; timepoint<=timepoint_start_WS_use+no_timesteps_WS_use; timepoint++) {
            // std::cout << "SAMPLING AT TIME: " << timepoint << std::endl;
            latt_traj->get_lattice_at_timepoint(timepoint)->wake_sleep_loop_rbm(i_opt_step, no_cd_steps, fname_coll.at(timepoint-timepoint_start_WS_use), options_wake_sleep);
        };
        
        clock_t t2 = clock();
        
        // Reap the moments
        // Before timepoint_start_A: don't slide means!
        // At & afer timepoint_start_A: slide means!
        // Never slide at timepoint_start_SIP (initial condition)
        
        for (auto timepoint=timepoint_start_WS_use; timepoint<=timepoint_start_WS_use+no_timesteps_WS_use; timepoint++) {
            latt_traj->get_lattice_at_timepoint(timepoint)->reap_ixn_moment_diffs();
        };
        
        // Print
        if (options.verbose) {
            for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
                
                // Print moment traj
                std::cout << ixn_param_traj->get_name() << " moments [" << timepoint_start_WS_use << "," << timepoint_start_WS_use+no_timesteps_WS_use << "]" << std::endl;
                ixn_param_traj->print_moment_diff_traj(timepoint_start_WS_use, no_timesteps_WS_use);
            };
        };
        
        /********************
         Solve diff eq for adjoint
         ********************/
        
        clock_t t3 = clock();
        
        // Set zero endpoint
        for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
            ixn_param_traj->get_adjoint()->set_timepoint_zero_end_cond(timepoint_start_A_use + no_timesteps_A_use);
        };
        
        if (options.l2_reg) {
            for (auto timepoint=timepoint_start_A_use + no_timesteps_A_use; timepoint>timepoint_start_A_use; timepoint--) {
                for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
                    if (!ixn_param_traj->get_is_val_fixed()) {
                        ixn_param_traj->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one_l2(timepoint,dt,options.l2_lambda.at(ixn_param_traj),options.l2_center.at(ixn_param_traj));
                    };
                };
            };
        } else if (options.l2_reg_traj) {
            for (auto timepoint=timepoint_start_A_use + no_timesteps_A_use; timepoint>timepoint_start_A_use; timepoint--) {
                for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
                    if (!ixn_param_traj->get_is_val_fixed()) {
                        ixn_param_traj->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one_l2(timepoint,dt,options.l2_lambda_traj.at(ixn_param_traj).at(timepoint),options.l2_center_traj.at(ixn_param_traj).at(timepoint));
                    };
                };
            };
        } else {
            for (auto timepoint=timepoint_start_A_use + no_timesteps_A_use; timepoint>timepoint_start_A_use; timepoint--) {
                for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
                    if (!ixn_param_traj->get_is_val_fixed()) {
                        ixn_param_traj->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one(timepoint,dt,form_abscissas);
                    };
                };
            };
        };
        
        /********************
         Form the update
         ********************/
        
        clock_t t4 = clock();
        
        for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
            if (!ixn_param_traj->get_is_val_fixed()) {
                ixn_param_traj->get_diff_eq_rhs()->update_calculate_and_store(timepoint_start_A_use,timepoint_start_A_use+no_timesteps_A_use,dt,form_abscissas);
            };
        };
        
        clock_t t5 = clock();
        
        if (options.verbose_timing) {
            double dt1 = (t1-t0)  / (double) CLOCKS_PER_SEC;
            double dt2 = (t2-t1)  / (double) CLOCKS_PER_SEC;
            double dt3 = (t3-t2)  / (double) CLOCKS_PER_SEC;
            double dt4 = (t4-t3)  / (double) CLOCKS_PER_SEC;
            double dt5 = (t5-t4)  / (double) CLOCKS_PER_SEC;
            double dt_tot = dt1 + dt2 + dt3 + dt4 + dt5;
            std::cout << "[time " << dt_tot << "] [F " << dt1/dt_tot << "] [wake/sleep " << dt2/dt_tot << "] [reap " << dt3/dt_tot << "] [adj " << dt4/dt_tot << "] [form update " << dt5/dt_tot << "]" << std::endl;
        };
    };
    
    void OptProblemDynamic::solve_one_step_rbm_params(std::shared_ptr<LatticeTrajAlternatingBinary> latt_traj, const std::vector<std::shared_ptr<IxnParamTraj>> &all_ixn_param_trajs, const std::vector<std::shared_ptr<Domain>> &all_domains, int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_cd_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_RBM options_wake_sleep) {
        
        solve_one_step_rbm_params_without_committ(latt_traj,all_ixn_param_trajs,all_domains,i_opt_step,timepoint_start_SIP,no_timesteps_SIP,timepoint_start_WS,no_timesteps_WS,timepoint_start_A,no_timesteps_A,dt,no_cd_steps,fname_traj_coll,options,options_wake_sleep);
        
        committ_step(latt_traj->get_all_ixn_param_trajs(), i_opt_step, options);
    };
    
    
    void OptProblemDynamic::solve_one_step_rbm_params_without_committ(std::shared_ptr<LatticeTrajAlternatingBinary> latt_traj, const std::vector<std::shared_ptr<IxnParamTraj>> &all_ixn_param_trajs, const std::vector<std::shared_ptr<Domain>> &all_domains, int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_cd_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_RBM options_wake_sleep) {
        
        if (latt_traj->get_no_markov_chains(MCType::AWAKE) != latt_traj->get_no_markov_chains(MCType::ASLEEP)) {
            std::cerr << ">>> OptProblemDynamic::solve_one_step_rbm_params_without_committ <<< asleep chains must equal awake chains" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        /*****
         Solve diff eq for F
         *****/
        
        clock_t t0 = clock();
        
        solve_ixn_param_trajs(all_ixn_param_trajs,all_domains, dt, timepoint_start_SIP, no_timesteps_SIP, options.no_steps_per_step_IP);
        
        // After now, no more need to reform coordinates
        bool form_abscissas = false;

        /*****
         Adjust adjoint range
         *****/
        
        int timepoint_start_A_use = timepoint_start_A;
        int no_timesteps_A_use = no_timesteps_A;
        int timepoint_start_WS_use = timepoint_start_WS;
        int no_timesteps_WS_use = no_timesteps_WS;
        
        if (options.locking_mode) {
            
            // Look at the cell in timepoint_start_A_use
            // Extend adjoint range to any timesteps also in this cell
            q3c1::Cell* cell = nullptr, *cell_prev = nullptr;
            bool repeat = false;
            if (timepoint_start_A_use != 0) {
                do {
                    // Dont repeat
                    repeat = false;
                    // Run through all ixn params
                    for (auto ixn: latt_traj->get_all_ixn_param_trajs()) {
                        // Get the cell at timepoint_start_A_use
                        cell = ixn->get_diff_eq_rhs()->get_cell_at_timepoint(timepoint_start_A_use);
                        cell_prev = ixn->get_diff_eq_rhs()->get_cell_at_timepoint(timepoint_start_A_use-1);
                        if (cell == cell_prev) {
                            // Same cell, extend adjoint domain
                            timepoint_start_A_use -= 1;
                            no_timesteps_A_use += 1;
                            // Ensure that this is covered by the wake-sleep range
                            if (timepoint_start_A_use < timepoint_start_WS_use) {
                                timepoint_start_WS_use = timepoint_start_A_use;
                            };
                            if (timepoint_start_WS_use + no_timesteps_WS_use < timepoint_start_A_use + no_timesteps_A_use) {
                                no_timesteps_WS_use = timepoint_start_A_use + no_timesteps_A_use - timepoint_start_WS_use;
                            };
                            // Restart the search
                            repeat = true;
                            break;
                        };
                    };
                } while (repeat && timepoint_start_A_use != 0);
            };
            std::cout << "Adjoint range is: " << timepoint_start_A_use << " to: " << timepoint_start_A_use + no_timesteps_A_use << std::endl;
            std::cout << "Wake/sleep range is: " << timepoint_start_WS_use << " to: " << timepoint_start_WS_use + no_timesteps_WS_use << std::endl;
            
            // Free all cells at timepoint start and after
            // Note: order matters, do this first!
            for (auto timepoint=timepoint_start_A_use; timepoint<=timepoint_start_A_use+no_timesteps_A_use; timepoint++) {
                for (auto ixn: latt_traj->get_all_ixn_param_trajs()) {
                    ixn->get_diff_eq_rhs()->fix_all_verts_around_at_timepoint(timepoint, false);
                };
            };
            
            // Fix all cells before the adjoint timepoint start
            // Note: order matters, do this second!
            for (auto timepoint=0; timepoint<timepoint_start_A_use; timepoint++) {
                for (auto ixn: latt_traj->get_all_ixn_param_trajs()) {
                    ixn->get_diff_eq_rhs()->fix_all_verts_around_at_timepoint(timepoint, true);
                };
            };
        };
        
        /*****
         Wake/asleep loop
         *****/
        
        clock_t t1 = clock();
        
        int no_awake_chains = latt_traj->get_no_markov_chains(MCType::AWAKE);
        std::vector<std::vector<FName>> fname_coll = fname_traj_coll.get_random_subset_fnames(no_awake_chains, timepoint_start_WS_use, no_timesteps_WS_use);
        
        for (auto timepoint=timepoint_start_WS_use; timepoint<=timepoint_start_WS_use+no_timesteps_WS_use; timepoint++) {
            latt_traj->get_lattice_at_timepoint(timepoint)->wake_sleep_loop_rbm(i_opt_step, no_cd_steps, fname_coll.at(timepoint-timepoint_start_WS_use), options_wake_sleep);
        };
        
        clock_t t2 = clock();
        
        // Reap the moments
        // Before timepoint_start_A: don't slide means!
        // At & afer timepoint_start_A: slide means!
        // Never slide at timepoint_start_SIP (initial condition)
        
        // NORMAL
        for (auto timepoint=timepoint_start_WS_use; timepoint<=timepoint_start_WS_use+no_timesteps_WS_use; timepoint++) {
            latt_traj->get_lattice_at_timepoint(timepoint)->reap_ixn_moment_diffs();
        };
        
        // Print
        if (options.verbose) {
            for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
                
                // Print moment traj
                std::cout << ixn_param_traj->get_name() << " moments [" << timepoint_start_WS_use << "," << timepoint_start_WS_use+no_timesteps_WS_use << "]" << std::endl;
                ixn_param_traj->print_moment_diff_traj(timepoint_start_WS_use, no_timesteps_WS_use);
            };
        };
        
        /********************
         Solve diff eq for adjoint
         ********************/
        
        clock_t t3 = clock();
        
        // Set zero endpoint
        for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
            ixn_param_traj->get_adjoint()->set_timepoint_zero_end_cond(timepoint_start_A_use + no_timesteps_A_use);
        };
        
        if (options.l2_reg) {
            for (auto timepoint=timepoint_start_A_use + no_timesteps_A_use; timepoint>timepoint_start_A_use; timepoint--) {
                for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
                    if (!ixn_param_traj->get_is_val_fixed()) {
                        ixn_param_traj->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one_l2(timepoint,dt,options.l2_lambda.at(ixn_param_traj),options.l2_center.at(ixn_param_traj));
                    };
                };
            };
        } else if (options.l2_reg_traj) {
            for (auto timepoint=timepoint_start_A_use + no_timesteps_A_use; timepoint>timepoint_start_A_use; timepoint--) {
                for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
                    if (!ixn_param_traj->get_is_val_fixed()) {
                        ixn_param_traj->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one_l2(timepoint,dt,options.l2_lambda_traj.at(ixn_param_traj).at(timepoint),options.l2_center_traj.at(ixn_param_traj).at(timepoint));
                    };
                };
            };
        } else {
            for (auto timepoint=timepoint_start_A_use + no_timesteps_A_use; timepoint>timepoint_start_A_use; timepoint--) {
                for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
                    if (!ixn_param_traj->get_is_val_fixed()) {
                        ixn_param_traj->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one(timepoint,dt,form_abscissas);
                    };
                };
            };
        };
        
        /********************
         Form the update
         ********************/
        
        clock_t t4 = clock();
        
        for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
            if (!ixn_param_traj->get_is_val_fixed()) {
                ixn_param_traj->get_diff_eq_rhs()->update_calculate_and_store(timepoint_start_A_use,timepoint_start_A_use+no_timesteps_A_use,dt, form_abscissas);
            };
        };
        
        clock_t t5 = clock();
        
        if (options.verbose_timing) {
            double dt1 = (t1-t0)  / (double) CLOCKS_PER_SEC;
            double dt2 = (t2-t1)  / (double) CLOCKS_PER_SEC;
            double dt3 = (t3-t2)  / (double) CLOCKS_PER_SEC;
            double dt4 = (t4-t3)  / (double) CLOCKS_PER_SEC;
            double dt5 = (t5-t4)  / (double) CLOCKS_PER_SEC;
            double dt_tot = dt1 + dt2 + dt3 + dt4 + dt5;
            std::cout << "[time " << dt_tot << "] [F " << dt1/dt_tot << "] [wake/sleep " << dt2/dt_tot << "] [reap " << dt3/dt_tot << "] [adj " << dt4/dt_tot << "] [form update " << dt5/dt_tot << "]" << std::endl;
        };
    };

    void OptProblemDynamic::solve_one_step_rbm_params(std::shared_ptr<LatticeTrajCenteredHom> latt_traj, const std::vector<std::shared_ptr<IxnParamTraj>> &all_ixn_param_trajs, const std::vector<std::shared_ptr<Domain>> &all_domains, int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_cd_steps, FNameTrajColl &fname_traj_coll, const std::vector<std::shared_ptr<AdjointParamsCenteredHomDerivTerm>> &all_deriv_terms, OptionsSolveDynamic options, OptionsWakeSleep_RBM options_wake_sleep) {
        
        solve_one_step_rbm_params_without_committ(latt_traj,all_ixn_param_trajs,all_domains,i_opt_step,timepoint_start_SIP,no_timesteps_SIP,timepoint_start_WS,no_timesteps_WS,timepoint_start_A,no_timesteps_A,dt,no_cd_steps,fname_traj_coll,all_deriv_terms,options,options_wake_sleep);
        
        committ_step(latt_traj->get_all_ixn_param_trajs(), i_opt_step, options);
    };
    
    
    void OptProblemDynamic::solve_one_step_rbm_params_without_committ(std::shared_ptr<LatticeTrajCenteredHom> latt_traj, const std::vector<std::shared_ptr<IxnParamTraj>> &all_ixn_param_trajs, const std::vector<std::shared_ptr<Domain>> &all_domains, int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_cd_steps, FNameTrajColl &fname_traj_coll, const std::vector<std::shared_ptr<AdjointParamsCenteredHomDerivTerm>> &all_deriv_terms, OptionsSolveDynamic options, OptionsWakeSleep_RBM options_wake_sleep) {
        
        if (latt_traj->get_no_markov_chains(MCType::AWAKE) != latt_traj->get_no_markov_chains(MCType::ASLEEP)) {
            std::cerr << ">>> OptProblemDynamic::solve_one_step_rbm_params_without_committ <<< asleep chains must equal awake chains" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        if (options.locking_mode) {
            std::cerr << ">>> OptProblemDynamic::solve_one_step_rbm_params_without_committ <<< Locking mode not supported here" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        if (options.l2_reg || options.l2_reg_traj) {
            std::cerr << ">>> OptProblemDynamic::solve_one_step_rbm_params_without_committ <<< L2 reg mode not supported here" << std::endl;
            exit(EXIT_FAILURE);
        };

        /*****
         Solve diff eq for F
         *****/
        
        clock_t t0 = clock();
        
        if (options.no_steps_per_step_IP == 1) {
            solve_ixn_param_trajs(all_ixn_param_trajs, all_domains, dt, timepoint_start_SIP, no_timesteps_SIP);
        } else {
            solve_ixn_param_trajs(all_ixn_param_trajs, all_domains, dt, timepoint_start_SIP, no_timesteps_SIP, options.no_steps_per_step_IP);
        };
        
        // After now, no more need to reform coordinates
        bool form_abscissas = false;

        /*****
         Wake/asleep loop
         *****/
        
        clock_t t1 = clock();
        
        int no_awake_chains = latt_traj->get_no_markov_chains(MCType::AWAKE);
        std::vector<std::vector<FName>> fname_coll = fname_traj_coll.get_random_subset_fnames(no_awake_chains, timepoint_start_WS, no_timesteps_WS);
        
        for (auto timepoint=timepoint_start_WS; timepoint<=timepoint_start_WS+no_timesteps_WS; timepoint++) {
            // std::cout << "SAMPLING AT TIME: " << timepoint << std::endl;
            latt_traj->get_lattice_at_timepoint(timepoint)->wake_sleep_loop_rbm(i_opt_step, no_cd_steps, fname_coll.at(timepoint-timepoint_start_WS), options_wake_sleep);
        };
        
        clock_t t2 = clock();
        
        // Reap the moments
        bool calculate_offset = false;
        for (auto timepoint=timepoint_start_WS; timepoint<=timepoint_start_WS+no_timesteps_WS; timepoint++) {
            latt_traj->get_lattice_centered_hom_at_timepoint(timepoint)->reap_ixn_moment_diffs_and_slide_centers(options.sliding_factor,calculate_offset);
        };
        
        // Print
        if (options.verbose) {
            for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
                std::cout << ixn_param_traj->get_name() << " moments [" << timepoint_start_WS << "," << timepoint_start_WS+no_timesteps_WS << "]" << std::endl;
                ixn_param_traj->print_moment_diff_traj(timepoint_start_WS, no_timesteps_WS);
            };
        };
        
        /********************
         Solve diff eq for adjoint
         ********************/
        
        clock_t t3 = clock();
        
        // Set zero endpoint
        for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
            ixn_param_traj->get_adjoint()->set_timepoint_zero_end_cond(timepoint_start_A + no_timesteps_A);
        };
        
        for (auto timepoint=timepoint_start_A + no_timesteps_A; timepoint>timepoint_start_A; timepoint--) {
            // Calculate deriv terms
            for (auto deriv_term: all_deriv_terms) {
                deriv_term->calculate_val_at_timepoint(timepoint,form_abscissas);
            };
            
            // Solve diff eq
            for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
                if (!ixn_param_traj->get_is_val_fixed()) {
                    ixn_param_traj->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one(timepoint, dt, form_abscissas);
                };
            };
        };
        
        /********************
         Form the update
         ********************/
        
        clock_t t4 = clock();
        
        for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
            if (!ixn_param_traj->get_is_val_fixed()) {
                ixn_param_traj->get_diff_eq_rhs()->update_calculate_and_store(timepoint_start_A,timepoint_start_A+no_timesteps_A,dt,form_abscissas);
            };
        };
        
        clock_t t5 = clock();
        
        if (options.verbose_timing) {
            double dt1 = (t1-t0)  / (double) CLOCKS_PER_SEC;
            double dt2 = (t2-t1)  / (double) CLOCKS_PER_SEC;
            double dt3 = (t3-t2)  / (double) CLOCKS_PER_SEC;
            double dt4 = (t4-t3)  / (double) CLOCKS_PER_SEC;
            double dt5 = (t5-t4)  / (double) CLOCKS_PER_SEC;
            double dt_tot = dt1 + dt2 + dt3 + dt4 + dt5;
            std::cout << "[time " << dt_tot << "] [F " << dt1/dt_tot << "] [wake/sleep " << dt2/dt_tot << "] [reap " << dt3/dt_tot << "] [adj " << dt4/dt_tot << "] [form update " << dt5/dt_tot << "]" << std::endl;
        };
    };

    
    
    // ***************
    // MARK: - RBM CD obs
    // ***************

    void OptProblemDynamic::solve_one_step_rbm_obs(std::shared_ptr<LatticeTraj> latt_traj, const std::vector<std::shared_ptr<IxnParamTraj>> &all_ixn_param_trajs, const std::vector<std::shared_ptr<Domain>> &all_domains, int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_cd_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_RBM options_wake_sleep) {
        
        solve_one_step_rbm_obs_without_committ(latt_traj,all_ixn_param_trajs,all_domains,i_opt_step,timepoint_start_SIP,no_timesteps_SIP,timepoint_start_WS,no_timesteps_WS,timepoint_start_A,no_timesteps_A,dt,no_cd_steps,fname_traj_coll,options,options_wake_sleep);
        
        committ_step(latt_traj->get_all_ixn_param_trajs(), i_opt_step, options);
    };
    
    
    void OptProblemDynamic::solve_one_step_rbm_obs_without_committ(std::shared_ptr<LatticeTraj> latt_traj, const std::vector<std::shared_ptr<IxnParamTraj>> &all_ixn_param_trajs, const std::vector<std::shared_ptr<Domain>> &all_domains, int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_cd_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_RBM options_wake_sleep) {
        
        if (latt_traj->get_no_markov_chains(MCType::AWAKE) != latt_traj->get_no_markov_chains(MCType::ASLEEP)) {
            std::cerr << ">>> OptProblemDynamic::solve_one_step_rbm_obs_without_committ <<< asleep chains must equal awake chains" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        if (options.locking_mode) {
            std::cerr << ">>> OptProblemDynamic::solve_one_step_rbm_obs_without_committ <<< locking mode not supported" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        if (options.l2_reg) {
            std::cerr << ">>> OptProblemDynamic::solve_one_step_rbm_obs_without_committ <<< l2_reg mode not supported" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        if (options.l2_reg_traj) {
            std::cerr << ">>> OptProblemDynamic::solve_one_step_rbm_obs_without_committ <<< l2_reg_traj mode not supported" << std::endl;
            exit(EXIT_FAILURE);
        };

        /*****
         Data for sampling
         *****/
        
        int no_awake_chains = latt_traj->get_no_markov_chains(MCType::AWAKE);
        std::vector<std::vector<FName>> fname_coll = fname_traj_coll.get_random_subset_fnames(no_awake_chains, timepoint_start_WS, no_timesteps_WS);
        
        /*****
         Get unique terms
         *****/
        
        std::vector<std::shared_ptr<AdjointObsCommonTerm>> all_common_terms;
        std::vector<std::shared_ptr<AdjointMomentCovTerm>> all_cov_terms;

        for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
            std::shared_ptr<AdjointObs> adjoint_obs = ixn_param_traj->get_adjoint_obs();
            std::vector< std::pair<std::shared_ptr<AdjointObsCommonTerm>,std::shared_ptr<AdjointMomentCovTerm>>> terms = adjoint_obs->get_terms();
            for (auto term_pr: terms) {
                auto it1 = std::find(all_common_terms.begin(), all_common_terms.end(), term_pr.first);
                if (it1 == all_common_terms.end()) {
                    all_common_terms.push_back(term_pr.first);
                };
                auto it2 = std::find(all_cov_terms.begin(), all_cov_terms.end(), term_pr.second);
                if (it2 == all_cov_terms.end()) {
                    all_cov_terms.push_back(term_pr.second);
                };
            };
        };
        
        /*****
         Get unique domains
         *****/
        
        std::vector<Domain1DObs*> all_domain_obs;

        for (auto domain: all_domains) {
            for (auto domain_obs: domain->get_domain_obs()) {
                auto it = std::find(all_domain_obs.begin(), all_domain_obs.end(), domain_obs);
                if (it == all_domain_obs.end()) {
                    all_domain_obs.push_back(domain_obs);
                };
            };
        };
        
        /*****
        Go through all timepoints
         *****/
        
        clock_t t0 = clock();

        std::shared_ptr<Lattice> latt;
        Iptr cov_ixn;
        double term_1, term_2, cross_term, domain_val;
        for (auto timepoint=timepoint_start_SIP; timepoint<=timepoint_start_SIP+no_timesteps_SIP; timepoint++) {
            
            latt = latt_traj->get_lattice_at_timepoint(timepoint);
            
            /*****
             Solve diff eq for F
             *****/
            
            if (timepoint != timepoint_start_SIP+no_timesteps_SIP) { // except at the very end
                solve_ixn_param_trajs_step(all_ixn_param_trajs,all_domains, dt, timepoint);
            };
            
            /*****
             Wake/asleep loop
             *****/
            
            if (timepoint >= timepoint_start_WS && timepoint<=timepoint_start_WS+no_timesteps_WS) {
               
               // Sample
                latt->wake_sleep_loop_rbm(i_opt_step, no_cd_steps, fname_coll.at(timepoint-timepoint_start_WS), options_wake_sleep);
                
                // Reap moments
                latt->reap_ixn_moment_diffs();
                
                // Reap extra moments for the adjoint system
                for (auto cov_term: all_cov_terms) {
                    cov_ixn = cov_term->get_ixn_param_traj()->get_ixn_param_at_timepoint(timepoint);
                    term_1 = latt->reap_moment(MCType::ASLEEP, cov_ixn);
                    term_2 = latt->reap_moment(MCType::ASLEEP, cov_term->get_layer_domain(), cov_term->get_species_domain());
                    cross_term = latt->reap_moment_adjoint_obs_cov_cross_term(cov_ixn, cov_term->get_layer_domain(), cov_term->get_species_domain());
                    
                    // Set
                    cov_term->set_val_1_at_timepoint(timepoint, cross_term);
                    cov_term->set_val_2_at_timepoint(timepoint, term_1);
                    cov_term->set_val_3_at_timepoint(timepoint, term_2);
                    
                    // std::cout << "COV TERMS: timepoint: " << timepoint << " cov_term: " << cov_term->get_ixn_param_traj()->get_name() << " " << cov_term->get_layer_domain() << " " << cov_term->get_species_domain()->get_name() << ": " << cross_term << " - " << term_1 << " * " << term_2 << std::endl;
                };
                
                // Get the new moments for the diff eq RHS
                for (auto domain_obs: all_domain_obs) {
                    domain_val = latt->reap_moment(MCType::ASLEEP, domain_obs->get_layer(), domain_obs->get_species());
                    domain_obs->set_val_at_timepoint(timepoint, domain_val);
                    
                    // std::cout << "DOMAIN: timepoint: " << timepoint << " domain: " << domain_obs->get_layer() << " " << domain_obs->get_species()->get_name() << ": " << domain_val << std::endl;
                };
            };
        };
        
        // Print both
        if (options.verbose) {
            for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
                std::cout << ixn_param_traj->get_name() << " [" << timepoint_start_WS << "," << timepoint_start_WS+no_timesteps_WS << "]" << std::endl;
                ixn_param_traj->print_val_traj(timepoint_start_WS, no_timesteps_WS);
                ixn_param_traj->print_moment_diff_traj(timepoint_start_WS, no_timesteps_WS);
            };
        };

        // Print ixn traj
        /*
        if (options.verbose) {
            for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
                std::cout << ixn_param_traj->get_name() << " ixns [" << timepoint_start_SIP << "," << timepoint_start_SIP+no_timesteps_SIP << "]" << std::endl;
                ixn_param_traj->print_val_traj(timepoint_start_SIP, no_timesteps_SIP);
            };
        };
        */
        
        // Print moment traj
        /*
        if (options.verbose) {
            for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
                std::cout << ixn_param_traj->get_name() << " moments [" << timepoint_start_WS << "," << timepoint_start_WS+no_timesteps_WS << "]" << std::endl;
                ixn_param_traj->print_moment_diff_traj(timepoint_start_WS, no_timesteps_WS);
            };
        };
        */
        
        // After now, no more need to reform coordinates
        bool form_abscissas = false;
        
        clock_t t1 = clock();

        /*****
         Solve diff eq for adjoint
         *****/
        
        // Set zero endpoint
        for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
            ixn_param_traj->get_adjoint()->set_timepoint_zero_end_cond(timepoint_start_A + no_timesteps_A);
        };
        
        for (auto timepoint=timepoint_start_A + no_timesteps_A; timepoint>timepoint_start_A; timepoint--) {
            // Evaluate the common terms
            for (auto common_term: all_common_terms) {
                // std::cout << "timepoint: " << timepoint << " common term: " << common_term->get_layer() << " " << common_term->get_species()->get_name() << std::endl;
                common_term->calculate_val_at_timepoint(timepoint);
            };
            
            // Solve diff eq
            for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
                if (!ixn_param_traj->get_is_val_fixed()) {
                    ixn_param_traj->get_adjoint_obs()->solve_diff_eq_at_timepoint_to_minus_one(timepoint,dt,form_abscissas);
                };
            };
        };
        
        /********************
         Form the update
         ********************/
        
        clock_t t2 = clock();
        
        for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
            if (!ixn_param_traj->get_is_val_fixed()) {
                ixn_param_traj->get_diff_eq_rhs()->update_calculate_and_store(timepoint_start_A,timepoint_start_A+no_timesteps_A,dt,form_abscissas);
            };
        };
        
        clock_t t3 = clock();
        
        if (options.verbose_timing) {
            double dt1 = (t1-t0)  / (double) CLOCKS_PER_SEC;
            double dt2 = (t2-t1)  / (double) CLOCKS_PER_SEC;
            double dt3 = (t3-t2)  / (double) CLOCKS_PER_SEC;
            double dt_tot = dt1 + dt2 + dt3;
            std::cout << "[time " << dt_tot << "] [F " << dt1/dt_tot << "] [wake/sleep " << dt2/dt_tot << "] [form update " << dt3/dt_tot << "]" << std::endl;
        };
    };
    
    // ***************
    // MARK: - 1D fully visible lattice
    // ***************
    
    void OptProblemDynamic::solve_one_step_1d_fully_visible(std::shared_ptr<LatticeTraj1DFullyVisible> latt_traj, const std::vector<std::shared_ptr<IxnParamTraj>> &all_ixn_param_trajs, const std::vector<std::shared_ptr<Domain>> &all_domains, int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_cd_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_1DFV_CD options_wake_sleep) {
        
        solve_one_step_1d_fully_visible_without_committ(latt_traj, all_ixn_param_trajs, all_domains, i_opt_step, timepoint_start_SIP, no_timesteps_SIP, timepoint_start_WS, no_timesteps_WS, timepoint_start_A, no_timesteps_A, dt, no_cd_steps, fname_traj_coll, options, options_wake_sleep);
        
        committ_step(latt_traj->get_all_ixn_param_trajs(), i_opt_step, options);
        
    };
    
    void OptProblemDynamic::solve_one_step_1d_fully_visible_without_committ(std::shared_ptr<LatticeTraj1DFullyVisible> latt_traj, const std::vector<std::shared_ptr<IxnParamTraj>> &all_ixn_param_trajs, const std::vector<std::shared_ptr<Domain>> &all_domains, int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_cd_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_1DFV_CD options_wake_sleep) {
        
        if (latt_traj->get_no_markov_chains(MCType::AWAKE) != latt_traj->get_no_markov_chains(MCType::ASLEEP)) {
            std::cerr << ">>> OptProblemDynamic::solve_one_step_1d_fully_visible_without_committ <<< asleep chains must equal awake chains" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        /*****
         Solve diff eq for F
         *****/
        
        clock_t t0 = clock();
        
        solve_ixn_param_trajs(all_ixn_param_trajs, all_domains, dt, timepoint_start_SIP, no_timesteps_SIP, options.no_steps_per_step_IP);
        
        // After now, no more need to reform coordinates
        bool form_abscissas = false;
        
        /*****
         Adjust adjoint range
         *****/
        
        int timepoint_start_A_use = timepoint_start_A;
        int no_timesteps_A_use = no_timesteps_A;
        int timepoint_start_WS_use = timepoint_start_WS;
        int no_timesteps_WS_use = no_timesteps_WS;
        
        if (options.locking_mode) {
            
            // Look at the cell in timepoint_start_A_use
            // Extend adjoint range to any timesteps also in this cell
            q3c1::Cell* cell = nullptr, *cell_prev = nullptr;
            bool repeat = false;
            if (timepoint_start_A_use != 0) {
                do {
                    // Dont repeat
                    repeat = false;
                    // Run through all ixn params
                    for (auto ixn: latt_traj->get_all_ixn_param_trajs()) {
                        // Get the cell at timepoint_start_A_use
                        cell = ixn->get_diff_eq_rhs()->get_cell_at_timepoint(timepoint_start_A_use);
                        cell_prev = ixn->get_diff_eq_rhs()->get_cell_at_timepoint(timepoint_start_A_use-1);
                        if (cell == cell_prev) {
                            // Same cell, extend adjoint domain
                            timepoint_start_A_use -= 1;
                            no_timesteps_A_use += 1;
                            // Ensure that this is covered by the wake-sleep range
                            if (timepoint_start_A_use < timepoint_start_WS_use) {
                                timepoint_start_WS_use = timepoint_start_A_use;
                            };
                            if (timepoint_start_WS_use + no_timesteps_WS_use < timepoint_start_A_use + no_timesteps_A_use) {
                                no_timesteps_WS_use = timepoint_start_A_use + no_timesteps_A_use - timepoint_start_WS_use;
                            };
                            // Restart the search
                            repeat = true;
                            break;
                        };
                    };
                } while (repeat && timepoint_start_A_use != 0);
            };
            std::cout << "Adjoint range is: " << timepoint_start_A_use << " to: " << timepoint_start_A_use + no_timesteps_A_use << std::endl;
            std::cout << "Wake/sleep range is: " << timepoint_start_WS_use << " to: " << timepoint_start_WS_use + no_timesteps_WS_use << std::endl;
            
            // Free all cells at timepoint start and after
            // Note: order matters, do this first!
            for (auto timepoint=timepoint_start_A_use; timepoint<=timepoint_start_A_use+no_timesteps_A_use; timepoint++) {
                for (auto ixn: latt_traj->get_all_ixn_param_trajs()) {
                    ixn->get_diff_eq_rhs()->fix_all_verts_around_at_timepoint(timepoint, false);
                };
            };
            
            // Fix all cells before the adjoint timepoint start
            // Note: order matters, do this second!
            for (auto timepoint=0; timepoint<timepoint_start_A_use; timepoint++) {
                for (auto ixn: latt_traj->get_all_ixn_param_trajs()) {
                    ixn->get_diff_eq_rhs()->fix_all_verts_around_at_timepoint(timepoint, true);
                };
            };
        };
        
        /*****
         Wake/asleep loop
         *****/
        
        clock_t t1 = clock();
        
        int no_awake_chains = latt_traj->get_no_markov_chains(MCType::AWAKE);
        std::vector<std::vector<FName>> fname_coll = fname_traj_coll.get_random_subset_fnames(no_awake_chains, timepoint_start_WS_use, no_timesteps_WS_use);
        
        for (auto timepoint=timepoint_start_WS_use; timepoint<=timepoint_start_WS_use+no_timesteps_WS_use; timepoint++) {
            latt_traj->get_lattice_at_timepoint(timepoint)->wake_sleep_loop_cd(i_opt_step, no_cd_steps, fname_coll.at(timepoint-timepoint_start_WS_use), options_wake_sleep);
        };
        
        clock_t t2 = clock();
        
        // Reap the moments
        // Before timepoint_start_A: don't slide means!
        // At & afer timepoint_start_A: slide means!
        // Never slide at timepoint_start_SIP (initial condition)
        
        for (auto timepoint=timepoint_start_WS_use; timepoint<=timepoint_start_WS_use+no_timesteps_WS_use; timepoint++) {
            latt_traj->get_lattice_at_timepoint(timepoint)->reap_ixn_moment_diffs();
        };
        
        // Print
        if (options.verbose) {
            for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
                
                // Print moment traj
                if (ixn_param_traj->get_type() == IxnParamType::W || ixn_param_traj->get_type() == IxnParamType::X) {
                    std::cout << ixn_param_traj->get_name() << " moments [" << timepoint_start_WS_use << "," << timepoint_start_WS_use+no_timesteps_WS_use << "]" << std::endl;
                    ixn_param_traj->print_moment_diff_traj(timepoint_start_WS_use, no_timesteps_WS_use);
                };
            };
        };
        
        /********************
         Solve diff eq for adjoint
         ********************/
        
        clock_t t3 = clock();
        
        // Set zero endpoint
        for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
            ixn_param_traj->get_adjoint()->set_timepoint_zero_end_cond(timepoint_start_A_use + no_timesteps_A_use);
        };
        
        if (options.l2_reg) {
            for (auto timepoint=timepoint_start_A_use + no_timesteps_A_use; timepoint>timepoint_start_A_use; timepoint--) {
                for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
                    if (!ixn_param_traj->get_is_val_fixed()) {
                        ixn_param_traj->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one_l2(timepoint,dt,options.l2_lambda.at(ixn_param_traj),options.l2_center.at(ixn_param_traj));
                    };
                };
            };
        } else {
            for (auto timepoint=timepoint_start_A_use + no_timesteps_A_use; timepoint>timepoint_start_A_use; timepoint--) {
                for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
                    if (!ixn_param_traj->get_is_val_fixed()) {
                        ixn_param_traj->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one(timepoint,dt,form_abscissas);
                    };
                };
            };
        };
        
        /********************
         Form the update
         ********************/
        
        clock_t t4 = clock();
        
        for (auto ixn_param_traj: latt_traj->get_all_ixn_param_trajs()) {
            if (!ixn_param_traj->get_is_val_fixed()) {
                ixn_param_traj->get_diff_eq_rhs()->update_calculate_and_store(timepoint_start_A_use,timepoint_start_A_use+no_timesteps_A_use,dt,form_abscissas);
            };
        };
        
        clock_t t5 = clock();
        
        if (options.verbose_timing) {
            double dt1 = (t1-t0)  / (double) CLOCKS_PER_SEC;
            double dt2 = (t2-t1)  / (double) CLOCKS_PER_SEC;
            double dt3 = (t3-t2)  / (double) CLOCKS_PER_SEC;
            double dt4 = (t4-t3)  / (double) CLOCKS_PER_SEC;
            double dt5 = (t5-t4)  / (double) CLOCKS_PER_SEC;
            double dt_tot = dt1 + dt2 + dt3 + dt4 + dt5;
            std::cout << "[time " << dt_tot << "] [F " << dt1/dt_tot << "] [wake/sleep " << dt2/dt_tot << "] [reap " << dt3/dt_tot << "] [adj " << dt4/dt_tot << "] [form update " << dt5/dt_tot << "]" << std::endl;
        };

        
    };
    
    void OptProblemDynamic::committ_step(const std::vector<std::shared_ptr<IxnParamTraj>> &ixn_params, int i_opt_step, OptionsSolveDynamic options) {
        
        if (options.solver == Solver::ADAM) {
            for (auto &ixn_param_traj: ixn_params) {
                if (!ixn_param_traj->get_is_val_fixed()) {
                    ixn_param_traj->get_diff_eq_rhs()->update_committ_stored_adam(i_opt_step,options.adam_beta_1,options.adam_beta_2,options.adam_eps);
                };
            };
        } else if (options.solver == Solver::SGD) {
            for (auto &ixn_param_traj: ixn_params) {
                if (!ixn_param_traj->get_is_val_fixed()) {
                    ixn_param_traj->get_diff_eq_rhs()->update_committ_stored_sgd();
                };
            };
        } else {
            std::cerr << ">>> OptProblemDynamic::solve_one_step <<< Solvers other than Adam are currently not supported" << std::endl;
            exit(EXIT_FAILURE);
        };
    };
};
