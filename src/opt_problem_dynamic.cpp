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
     Init structures
     ********************/

    /*
    void OptProblemDynamic::set_no_markov_chains(MCType chain, int no_markov_chains) {
        
        // Store
        _no_markov_chains[chain] = no_markov_chains;
     
        // Lattice
        _latt_traj->set_no_markov_chains(chain, no_markov_chains);
    };
    void OptProblemDynamic::set_no_timesteps_ixn_params(int no_timesteps_ixn_params) {
        _no_timesteps_ixn_params = no_timesteps_ixn_params;
        _no_timepoints_ixn_params = _no_timesteps_ixn_params + 1;
        for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
            ixn_param_traj->set_no_timesteps(_no_timesteps_ixn_params);
        };
    };
    void OptProblemDynamic::set_no_timesteps_lattice(int timepoint_start_lattice, int no_timesteps_lattice) {
        _timepoint_start_lattice = timepoint_start_lattice;
        _no_timesteps_lattice = no_timesteps_lattice;
        _no_timepoints_lattice = _no_timesteps_lattice + 1;
        
        _latt_traj->set_no_timesteps(_timepoint_start_lattice, _no_timesteps_lattice);
    };
    */
    
    /********************
     Solve
     ********************/
    
    // Check if options passed are valid
    void OptProblemDynamic::check_options(int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_mean_field_updates, int no_gibbs_sampling_steps, OptionsSolveDynamic options, OptionsWakeSleep options_wake_sleep) {
        if (options.solver == Solver::SGD || options.solver == Solver::NESTEROV) {
            std::cerr << ">>> OptProblemDynamic::check_options <<< Error: only Adam is currently supported as a solver" << std::endl;
            exit(EXIT_FAILURE);
        };
    };
    
    // One step
    void OptProblemDynamic::solve_one_step(int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_mean_field_updates, int no_gibbs_sampling_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep options_wake_sleep) {
        
        /*****
         Check options
         *****/
        
        if (options.should_check_options) {
            check_options(timepoint_start_SIP, no_timesteps_SIP, timepoint_start_WS, no_timesteps_WS, timepoint_start_A, no_timesteps_A, dt, no_mean_field_updates, no_gibbs_sampling_steps, options, options_wake_sleep);
        };
        
        /*****
         Solve diff eq for F
         *****/
        
        // Solve over all time
        std::map<ITptr,double> vals;
        std::map<ITptr,double> vals_next;
        for (auto timepoint=timepoint_start_SIP; timepoint<timepoint_start_SIP+no_timesteps_SIP; timepoint++) {
            // Get initial vals at this timepoint
            for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                vals[ixn_param_traj] = ixn_param_traj->get_ixn_param_at_timepoint(timepoint)->get_val();
            };
            
            // Calculate diff eqs to a new step
            for (auto i=0; i<options.no_steps_per_step_IP; i++) {
                for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                    vals_next[ixn_param_traj] = vals.at(ixn_param_traj) + (dt / options.no_steps_per_step_IP) * ixn_param_traj->get_diff_eq_rhs()->get_val_from_map(vals);
                };
                
                // Advance
                vals = vals_next;
            };
            
            // Write final vals
            for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                ixn_param_traj->get_ixn_param_at_timepoint(timepoint+1)->set_val(vals.at(ixn_param_traj));
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
        
        /*****
         Wake/asleep loop
         *****/
        
        int no_awake_chains = _latt_traj->get_no_markov_chains(MCType::AWAKE);
        std::vector<std::vector<FName>> fname_coll = fname_traj_coll.get_random_subset_fnames(no_awake_chains, timepoint_start_WS, no_timesteps_WS);
        
        for (auto timepoint=timepoint_start_WS; timepoint<=timepoint_start_WS+no_timesteps_WS; timepoint++) {
            
            _latt_traj->get_lattice_at_timepoint(timepoint)->wake_sleep_loop(i_opt_step, no_mean_field_updates, no_gibbs_sampling_steps, fname_coll.at(timepoint-timepoint_start_WS), options_wake_sleep);
            
        };
        
        // Print
        if (options.verbose) {
            for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                
                // Print traj of ixn params
                std::cout << ixn_param_traj->get_name() << " [" << timepoint_start_SIP << "," << timepoint_start_SIP+no_timesteps_SIP << "]" << std::endl;
                ixn_param_traj->print_val_traj(timepoint_start_SIP, no_timesteps_SIP);
                
                // Print moment traj
                std::cout << ixn_param_traj->get_name() << " moments [" << timepoint_start_WS << "," << timepoint_start_WS+no_timesteps_WS << "]" << std::endl;
                ixn_param_traj->print_moment_traj(timepoint_start_WS, no_timesteps_WS);
            };
        };
        
        /********************
         Solve diff eq for adjoint
         ********************/
        
        // Set zero endpoint
        for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
            ixn_param_traj->get_adjoint()->set_timepoint_zero_end_cond(timepoint_start_A + no_timesteps_A);
        };
        
        if (options.l2_reg) {
            for (auto timepoint=timepoint_start_A + no_timesteps_A; timepoint>timepoint_start_A; timepoint--) {
                for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                    if (!ixn_param_traj->get_is_val_fixed_to_init_cond() && !ixn_param_traj->get_are_vals_fixed()) {
                        ixn_param_traj->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one_l2(timepoint,dt,options.l2_lambda,options.l2_center);
                    };
                };
            };
        } else {
            for (auto timepoint=timepoint_start_A + no_timesteps_A; timepoint>timepoint_start_A; timepoint--) {
                for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                    if (!ixn_param_traj->get_is_val_fixed_to_init_cond() && !ixn_param_traj->get_are_vals_fixed()) {
                        ixn_param_traj->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one(timepoint,dt);
                    };
                };
            };
        };

        /********************
         Form the update
         ********************/

        for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
            if (!ixn_param_traj->get_is_val_fixed_to_init_cond() && !ixn_param_traj->get_are_vals_fixed()) {
                
                ixn_param_traj->get_diff_eq_rhs()->update_calculate_and_store(timepoint_start_A,timepoint_start_A+no_timesteps_A,dt);
            };
        };
        
        /********************
         Committ the update
         ********************/
        
        if (options.solver == Solver::ADAM) {
            for (auto &ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                if (!ixn_param_traj->get_is_val_fixed_to_init_cond() && !ixn_param_traj->get_are_vals_fixed()) {
                    ixn_param_traj->get_diff_eq_rhs()->update_committ_stored_adam(i_opt_step,options.adam_beta_1,options.adam_beta_2,options.adam_eps);
                };
            };
        } else {
            std::cerr << ">>> OptProblemDynamic::solve_one_step <<< Solvers other than Adam are currently not supported" << std::endl;
            exit(EXIT_FAILURE);
        };
    };
};
