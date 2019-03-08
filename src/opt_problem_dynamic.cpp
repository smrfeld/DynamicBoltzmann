#include "../include/dblz_bits/opt_problem_dynamic.hpp"

// Other headers
#include "../include/dblz_bits/ixn_param.hpp"
#include "../include/dblz_bits/ixn_param_traj.hpp"
#include "../include/dblz_bits/lattice.hpp"
#include "../include/dblz_bits/lattice_traj.hpp"
#include "../include/dblz_bits/moment.hpp"
#include "../include/dblz_bits/fname_traj.hpp"
#include "../include/dblz_bits/general.hpp"
#include "../include/dblz_bits/adjoint.hpp"
#include "../include/dblz_bits/diff_eq_rhs.hpp"

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
    
    OptProblemDynamic::OptProblemDynamic(std::shared_ptr<LatticeTraj> latt_traj) {
        _latt_traj = latt_traj;
    };
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
    void OptProblemDynamic::_move(OptProblemDynamic &other) {
        _latt_traj = std::move(other._latt_traj);
    };
    void OptProblemDynamic::_copy(const OptProblemDynamic& other) {
        if (other._latt_traj) {
            _latt_traj = std::make_shared<LatticeTraj>(*other._latt_traj);
        } else {
            _latt_traj = nullptr;
        };
    };
    
    /********************
     Solve
     ********************/
    
    // Solve helpers
    void OptProblemDynamic::solve_ixn_param_trajs(double dt, int timepoint_start, int no_timesteps, int no_steps_per_step) const {
        
        // Solve over all time
        std::map<ITptr,std::vector<double>> vals;
        for (auto timepoint=timepoint_start; timepoint<timepoint_start+no_timesteps; timepoint++) {
            // Get initial vals at this timepoint
            for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                vals[ixn_param_traj] = std::vector<double>(no_steps_per_step+1);
                vals[ixn_param_traj][0] = ixn_param_traj->get_ixn_param_at_timepoint(timepoint)->get_val();
            };
            
            // Calculate diff eqs to a new step
            for (auto i=0; i<no_steps_per_step; i++) {
                for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                    if (!ixn_param_traj->get_is_val_fixed()) {
                        // std::cout << ixn_param_traj->get_name() << " " << timepoint << " " << i << " " << ixn_param_traj->get_diff_eq_rhs()->get_val_from_map(vals, i) << std::endl;
                        
                        vals[ixn_param_traj][i+1] = vals.at(ixn_param_traj).at(i) + (dt / no_steps_per_step) * ixn_param_traj->get_diff_eq_rhs()->get_val_from_map(vals, i);
                    } else {
                        vals[ixn_param_traj][i+1] = vals.at(ixn_param_traj).at(i);
                    };
                };
            };
            
            // Write final vals
            for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                if (!ixn_param_traj->get_is_val_fixed()) {
                    ixn_param_traj->get_ixn_param_at_timepoint(timepoint+1)->set_val(vals.at(ixn_param_traj).at(no_steps_per_step));
                };
            };
        };
        
        /*
         for (auto timepoint=timepoint_start_SIP; timepoint<timepoint_start_SIP+no_timesteps_SIP; timepoint++) {
         for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
         if (!ixn_param_traj->get_is_val_fixed_to_init_cond() && !ixn_param_traj->get_are_vals_fixed()) {
         ixn_param_traj->solve_diff_eq_at_timepoint_to_plus_one(timepoint,dt);
         };
         };
         };
         */
    };

    
    // One step
    void OptProblemDynamic::solve_one_step_bm_pcd(int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_mean_field_updates, int no_gibbs_sampling_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_BM_PCD options_wake_sleep) {

    solve_one_step_bm_pcd_without_committ(i_opt_step,timepoint_start_SIP,no_timesteps_SIP,timepoint_start_WS,no_timesteps_WS,timepoint_start_A,no_timesteps_A,dt,no_mean_field_updates,no_gibbs_sampling_steps,fname_traj_coll,options,options_wake_sleep);
        
        committ_step(i_opt_step, options);
    };
    
    
    void OptProblemDynamic::solve_one_step_bm_pcd_without_committ(int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_mean_field_updates, int no_gibbs_sampling_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_BM_PCD options_wake_sleep) {
        
        /*****
         Solve diff eq for F
         *****/
        
        clock_t t0 = clock();
        
        solve_ixn_param_trajs(dt, timepoint_start_SIP, no_timesteps_SIP, options.no_steps_per_step_IP);
        
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
                    for (auto ixn: _latt_traj->get_all_ixn_param_trajs()) {
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
            
            // Fix all cells before the adjoint timepoint start
            for (auto timepoint=0; timepoint<timepoint_start_A_use; timepoint++) {
                for (auto ixn: _latt_traj->get_all_ixn_param_trajs()) {
                    ixn->get_diff_eq_rhs()->fix_all_verts_around_at_timepoint(timepoint, true);
                };
            };
            
            // Free all cells at timepoint start and after
            for (auto timepoint=timepoint_start_A_use; timepoint<=timepoint_start_A_use+no_timesteps_A_use; timepoint++) {
                for (auto ixn: _latt_traj->get_all_ixn_param_trajs()) {
                    ixn->get_diff_eq_rhs()->fix_all_verts_around_at_timepoint(timepoint, false);
                };
            };
        };
        
        /*****
         Wake/asleep loop
         *****/
        
        clock_t t1 = clock();
        
        int no_awake_chains = _latt_traj->get_no_markov_chains(MCType::AWAKE);
        std::vector<std::vector<FName>> fname_coll = fname_traj_coll.get_random_subset_fnames(no_awake_chains, timepoint_start_WS_use, no_timesteps_WS_use);

        for (auto timepoint=timepoint_start_WS_use; timepoint<=timepoint_start_WS_use+no_timesteps_WS_use; timepoint++) {
            _latt_traj->get_lattice_at_timepoint(timepoint)->wake_sleep_loop_bm_pcd(i_opt_step, no_mean_field_updates, no_gibbs_sampling_steps, fname_coll.at(timepoint-timepoint_start_WS_use), options_wake_sleep);
        };
        
        clock_t t2 = clock();
        
        // Reap the moments
        // Before timepoint_start_A: don't slide means!
        // At & afer timepoint_start_A: slide means!
        // Never slide at timepoint_start_SIP (initial condition)
        bool slide_means;
        for (auto timepoint=timepoint_start_WS_use; timepoint<=timepoint_start_WS_use+no_timesteps_WS_use; timepoint++) {
            if (timepoint < timepoint_start_A_use) {
                slide_means = false;
            } else {
                slide_means = true;
            };
            if (timepoint == timepoint_start_SIP) {
                slide_means = false;
            };
            
            if (_latt_traj->get_lattice_at_timepoint(timepoint)->get_lattice_mode() == LatticeMode::NORMAL) {
                _latt_traj->get_lattice_at_timepoint(timepoint)->reap_moments_and_slide_centers_normal();
            } else if (_latt_traj->get_lattice_at_timepoint(timepoint)->get_lattice_mode() == LatticeMode::NORMAL_W_CENTERED_GRADIENT_PT) {
                _latt_traj->get_lattice_at_timepoint(timepoint)->reap_moments_and_slide_centers_normal_w_centered_gradient_pt(slide_means,options.sliding_factor);
            } else if (_latt_traj->get_lattice_at_timepoint(timepoint)->get_lattice_mode() == LatticeMode::NORMAL_W_CENTERED_GRADIENT_VEC) {
                _latt_traj->get_lattice_at_timepoint(timepoint)->reap_moments_and_slide_centers_normal_w_centered_gradient_vec(slide_means,options.sliding_factor);
            } else if (_latt_traj->get_lattice_at_timepoint(timepoint)->get_lattice_mode() == LatticeMode::CENTERED_PT) {
                _latt_traj->get_lattice_at_timepoint(timepoint)->reap_moments_and_slide_centers_centered_pt(slide_means,options.sliding_factor);
            };
        };
        
        // Print
        if (options.verbose) {
            for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                
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
        for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
            ixn_param_traj->get_adjoint()->set_timepoint_zero_end_cond(timepoint_start_A_use + no_timesteps_A_use);
        };
        
        if (options.l2_reg) {
            if (options.l2_reg_center_traj) {
                for (auto timepoint=timepoint_start_A_use + no_timesteps_A_use; timepoint>timepoint_start_A_use; timepoint--) {
                    for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                        if (!ixn_param_traj->get_is_val_fixed()) {
                            ixn_param_traj->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one_l2(timepoint,dt,options.l2_lambda.at(ixn_param_traj),options.l2_center_traj.at(ixn_param_traj).at(timepoint));
                        };
                    };
                };
            } else {
                for (auto timepoint=timepoint_start_A_use + no_timesteps_A_use; timepoint>timepoint_start_A_use; timepoint--) {
                    for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                        if (!ixn_param_traj->get_is_val_fixed()) {
                            ixn_param_traj->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one_l2(timepoint,dt,options.l2_lambda.at(ixn_param_traj),options.l2_center.at(ixn_param_traj));
                        };
                    };
                };
            };
        } else {
            for (auto timepoint=timepoint_start_A_use + no_timesteps_A_use; timepoint>timepoint_start_A_use; timepoint--) {
                for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                    if (!ixn_param_traj->get_is_val_fixed()) {
                        ixn_param_traj->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one(timepoint,dt);
                    };
                };
            };
        };
        
        /********************
         Form the update
         ********************/
        
        clock_t t4 = clock();
        
        for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
            if (!ixn_param_traj->get_is_val_fixed()) {
                ixn_param_traj->get_diff_eq_rhs()->update_calculate_and_store(timepoint_start_A_use,timepoint_start_A_use+no_timesteps_A_use,dt);
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
    
    void OptProblemDynamic::solve_one_step_rbm_cd(int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_cd_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_RBM_CD options_wake_sleep) {
        
        solve_one_step_rbm_cd_without_committ(i_opt_step,timepoint_start_SIP,no_timesteps_SIP,timepoint_start_WS,no_timesteps_WS,timepoint_start_A,no_timesteps_A,dt,no_cd_steps,fname_traj_coll,options,options_wake_sleep);
        
        committ_step(i_opt_step, options);
    };
    
    
    void OptProblemDynamic::solve_one_step_rbm_cd_without_committ(int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_cd_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_RBM_CD options_wake_sleep) {
        
        if (_latt_traj->get_no_markov_chains(MCType::AWAKE) != _latt_traj->get_no_markov_chains(MCType::ASLEEP)) {
            std::cerr << ">>> OptProblemDynamic::solve_one_step_rbm_cd_without_committ <<< asleep chains must equal awake chains" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        /*****
         Solve diff eq for F
         *****/
        
        clock_t t0 = clock();
        
        solve_ixn_param_trajs(dt, timepoint_start_SIP, no_timesteps_SIP, options.no_steps_per_step_IP);

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
                    for (auto ixn: _latt_traj->get_all_ixn_param_trajs()) {
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
            
            // Fix all cells before the adjoint timepoint start
            for (auto timepoint=0; timepoint<timepoint_start_A_use; timepoint++) {
                for (auto ixn: _latt_traj->get_all_ixn_param_trajs()) {
                    ixn->get_diff_eq_rhs()->fix_all_verts_around_at_timepoint(timepoint, true);
                };
            };
            
            // Free all cells at timepoint start and after
            for (auto timepoint=timepoint_start_A_use; timepoint<=timepoint_start_A_use+no_timesteps_A_use; timepoint++) {
                for (auto ixn: _latt_traj->get_all_ixn_param_trajs()) {
                    ixn->get_diff_eq_rhs()->fix_all_verts_around_at_timepoint(timepoint, false);
                };
            };
        };
        
        /*****
         Wake/asleep loop
         *****/
        
        clock_t t1 = clock();
        
        int no_awake_chains = _latt_traj->get_no_markov_chains(MCType::AWAKE);
        std::vector<std::vector<FName>> fname_coll = fname_traj_coll.get_random_subset_fnames(no_awake_chains, timepoint_start_WS_use, no_timesteps_WS_use);
        
        for (auto timepoint=timepoint_start_WS_use; timepoint<=timepoint_start_WS_use+no_timesteps_WS_use; timepoint++) {
            _latt_traj->get_lattice_at_timepoint(timepoint)->wake_sleep_loop_rbm_cd(i_opt_step, no_cd_steps, fname_coll.at(timepoint-timepoint_start_WS_use), options_wake_sleep);
        };
        
        clock_t t2 = clock();
        
        // Reap the moments
        // Before timepoint_start_A: don't slide means!
        // At & afer timepoint_start_A: slide means!
        // Never slide at timepoint_start_SIP (initial condition)
        
        if (_latt_traj->get_lattice_mode() == LatticeMode::NORMAL) {
            // NORMAL
            for (auto timepoint=timepoint_start_WS_use; timepoint<=timepoint_start_WS_use+no_timesteps_WS_use; timepoint++) {
                _latt_traj->get_lattice_at_timepoint(timepoint)->reap_moments_and_slide_centers_normal();
            };
        } else if (_latt_traj->get_lattice_mode() == LatticeMode::NORMAL_W_CENTERED_GRADIENT_PT) {
            // NORMAL W CENTERED GRADIENT PT
            bool slide_means;
            for (auto timepoint=timepoint_start_WS_use; timepoint<=timepoint_start_WS_use+no_timesteps_WS_use; timepoint++) {
                if (timepoint < timepoint_start_A_use) {
                    slide_means = false;
                } else {
                    slide_means = true;
                };
                if (timepoint == timepoint_start_SIP) {
                    slide_means = false;
                };
                _latt_traj->get_lattice_at_timepoint(timepoint)->reap_moments_and_slide_centers_normal_w_centered_gradient_pt(slide_means,options.sliding_factor);
            };
        } else if (_latt_traj->get_lattice_mode() == LatticeMode::NORMAL_W_CENTERED_GRADIENT_VEC) {
            // NORMAL W CENTERED GRADIENT VEC
            bool slide_means;
            for (auto timepoint=timepoint_start_WS_use; timepoint<=timepoint_start_WS_use+no_timesteps_WS_use; timepoint++) {
                if (timepoint < timepoint_start_A_use) {
                    slide_means = false;
                } else {
                    slide_means = true;
                };
                if (timepoint == timepoint_start_SIP) {
                    slide_means = false;
                };
                _latt_traj->get_lattice_at_timepoint(timepoint)->reap_moments_and_slide_centers_normal_w_centered_gradient_vec(slide_means,options.sliding_factor);
            };
        } else if (_latt_traj->get_lattice_mode() == LatticeMode::CENTERED_PT) {
            // CENTERED PT
            bool slide_means;
            for (auto timepoint=timepoint_start_WS_use; timepoint<=timepoint_start_WS_use+no_timesteps_WS_use; timepoint++) {
                if (timepoint < timepoint_start_A_use) {
                    slide_means = false;
                } else {
                    slide_means = true;
                };
                if (timepoint == timepoint_start_SIP) {
                    slide_means = false;
                };
                _latt_traj->get_lattice_at_timepoint(timepoint)->reap_moments_and_slide_centers_centered_pt(slide_means,options.sliding_factor);
            };
        };
        
        // Print
        if (options.verbose) {
            for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                
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
        for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
            ixn_param_traj->get_adjoint()->set_timepoint_zero_end_cond(timepoint_start_A_use + no_timesteps_A_use);
        };
        
        if (options.l2_reg) {
            if (options.l2_reg_center_traj) {
                for (auto timepoint=timepoint_start_A_use + no_timesteps_A_use; timepoint>timepoint_start_A_use; timepoint--) {
                    for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                        if (!ixn_param_traj->get_is_val_fixed()) {
                            ixn_param_traj->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one_l2(timepoint,dt,options.l2_lambda.at(ixn_param_traj),options.l2_center_traj.at(ixn_param_traj).at(timepoint));
                        };
                    };
                };
            } else {
                for (auto timepoint=timepoint_start_A_use + no_timesteps_A_use; timepoint>timepoint_start_A_use; timepoint--) {
                    for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                        if (!ixn_param_traj->get_is_val_fixed()) {
                            ixn_param_traj->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one_l2(timepoint,dt,options.l2_lambda.at(ixn_param_traj),options.l2_center.at(ixn_param_traj));
                        };
                    };
                };
            };
        } else {
            for (auto timepoint=timepoint_start_A_use + no_timesteps_A_use; timepoint>timepoint_start_A_use; timepoint--) {
                for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                    if (!ixn_param_traj->get_is_val_fixed()) {
                        ixn_param_traj->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one(timepoint,dt);
                    };
                };
            };
        };
        
        /********************
         Form the update
         ********************/
        
        clock_t t4 = clock();
        
        for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
            if (!ixn_param_traj->get_is_val_fixed()) {
                ixn_param_traj->get_diff_eq_rhs()->update_calculate_and_store(timepoint_start_A_use,timepoint_start_A_use+no_timesteps_A_use,dt);
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

    
    void OptProblemDynamic::committ_step(int i_opt_step, OptionsSolveDynamic options) {
        
        if (options.solver == Solver::ADAM) {
            for (auto &ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                if (!ixn_param_traj->get_is_val_fixed()) {
                    ixn_param_traj->get_diff_eq_rhs()->update_committ_stored_adam(i_opt_step,options.adam_beta_1,options.adam_beta_2,options.adam_eps);
                };
            };
        } else if (options.solver == Solver::SGD) {
            for (auto &ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
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
