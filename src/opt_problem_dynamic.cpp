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
    
    OptProblemDynamic::OptProblemDynamic(std::shared_ptr<LatticeTraj> latt_traj, int no_markov_chains_awake, int no_markov_chains_asleep, int no_timesteps_ixn_params, int timepoint_start_lattice, int no_timesteps_lattice) {
        _latt_traj = latt_traj;

        set_no_timesteps_ixn_params(no_timesteps_ixn_params); // NOTE: must set this first!
        set_no_timesteps_lattice(timepoint_start_lattice, no_timesteps_lattice);
        set_no_markov_chains(MCType::AWAKE,no_markov_chains_awake);
        set_no_markov_chains(MCType::ASLEEP,no_markov_chains_asleep);
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
        _no_markov_chains = other._no_markov_chains;
        _latt_traj = std::move(other._latt_traj);
        _no_timesteps_ixn_params = other._no_timesteps_ixn_params;
        _no_timepoints_ixn_params = other._no_timepoints_ixn_params;
        _timepoint_start_lattice = other._timepoint_start_lattice;
        _no_timesteps_lattice = other._no_timesteps_lattice;
        _no_timepoints_lattice = other._no_timepoints_lattice;
    };
    void OptProblemDynamic::_copy(const OptProblemDynamic& other) {
        _no_markov_chains = other._no_markov_chains;
        if (other._latt_traj) {
            _latt_traj = std::make_shared<LatticeTraj>(*other._latt_traj);
        };
        _no_timesteps_ixn_params = other._no_timesteps_ixn_params;
        _no_timepoints_ixn_params = other._no_timepoints_ixn_params;
        _timepoint_start_lattice = other._timepoint_start_lattice;
        _no_timesteps_lattice = other._no_timesteps_lattice;
        _no_timepoints_lattice = other._no_timepoints_lattice;
    };
    
    /********************
     Init structures
     ********************/

    void OptProblemDynamic::set_no_markov_chains(MCType chain, int no_markov_chains) {
        
        // Store
        _no_markov_chains[chain] = no_markov_chains;
        
        // Moments
        /*
        for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
            for (auto timepoint=0; timepoint<_no_timepoints_ixn_params; timepoint++) {
                ixn_param_traj->get_ixn_param_at_timepoint(timepoint)->get_moment()->set_no_markov_chains(chain, no_markov_chains);
            };
        };
         */
        
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
    
    /********************
     Solve
     ********************/
    
    // Check if options passed are valid
    void OptProblemDynamic::check_options(int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WSA, int no_timesteps_WSA, double dt, int no_mean_field_updates, int no_gibbs_sampling_steps, OptionsSolveDynamic options, OptionsWakeSleep options_wake_sleep) {
        if (options.solver == Solver::SGD || options.solver == Solver::NESTEROV) {
            std::cerr << ">>> OptProblemDynamic::check_options <<< Error: only Adam is currently supported as a solver" << std::endl;
            exit(EXIT_FAILURE);
        };
    };
    
    // One step
    void OptProblemDynamic::solve_one_step(int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WSA, int no_timesteps_WSA, double dt, int no_mean_field_updates, int no_gibbs_sampling_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep options_wake_sleep) {
        
        /*****
         Check options
         *****/
        
        if (options.verbose) {
            std::cout << "--- Checking options ---" << std::endl;
        };
        
        if (options.should_check_options) {
            check_options(timepoint_start_SIP, no_timesteps_SIP, timepoint_start_WSA, no_timesteps_WSA, dt, no_mean_field_updates, no_gibbs_sampling_steps, options, options_wake_sleep);
        };
        
        if (options.verbose) {
            std::cout << "--- [Finished] ---" << std::endl;
        };
        
        /*****
         Solve diff eq for F
         *****/
        
        if (options.verbose) {
            std::cout << "--- Solving diff eqs ---" << std::endl;
        };

        // Solve over all time
        for (auto timepoint=timepoint_start_SIP; timepoint<timepoint_start_SIP+no_timesteps_SIP; timepoint++) {
            for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                if (!ixn_param_traj->get_is_val_fixed_to_init_cond() && !ixn_param_traj->get_are_vals_fixed()) {
                    ixn_param_traj->solve_diff_eq_at_timepoint_to_plus_one(timepoint,dt);
                };
            };
        };
        
        if (options.verbose) {
            std::cout << "--- [Finished] ---" << std::endl;
        };

        /*****
         Wake/asleep loop
         *****/
        
        if (options.verbose) {
            std::cout << "--- Wake-sleep ---" << std::endl;
        };
        
        std::vector<std::vector<FName>> fname_coll = fname_traj_coll.get_random_subset_fnames(_no_markov_chains.at(MCType::AWAKE), timepoint_start_WSA, no_timesteps_WSA);
        
        for (auto timepoint=timepoint_start_WSA; timepoint<=timepoint_start_WSA+no_timesteps_WSA; timepoint++) {
            
            _latt_traj->get_lattice_at_timepoint(timepoint)->wake_sleep_loop(i_opt_step, no_mean_field_updates, no_gibbs_sampling_steps, fname_coll.at(timepoint-timepoint_start_WSA), options_wake_sleep);
            
        };
        
        // Print
        for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
            std::cout << ixn_param_traj->get_name() << " [" << _timepoint_start_lattice << "," << _timepoint_start_lattice+_no_timesteps_lattice << "]" << std::endl;
            
            // Print traj of ixn params
            ixn_param_traj->print_val_traj(_timepoint_start_lattice, _no_timesteps_lattice);
            
            // Print moment traj
            ixn_param_traj->print_moment_traj(_timepoint_start_lattice, _no_timesteps_lattice);
        };
        
        if (options.verbose) {
            std::cout << "--- [Finished] ---" << std::endl;
        };

        /********************
         Solve diff eq for adjoint
         ********************/
        
        if (options.verbose) {
            std::cout << "--- Solving adjoints ---" << std::endl;
        };

        // Set zero endpoint
        for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
            ixn_param_traj->get_adjoint()->set_timepoint_zero_end_cond(timepoint_start_WSA + no_timesteps_WSA);
        };
        
        if (options.l2_reg) {
            for (auto timepoint=timepoint_start_WSA + no_timesteps_WSA; timepoint>timepoint_start_WSA; timepoint--) {
                for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                    if (!ixn_param_traj->get_is_val_fixed_to_init_cond() && !ixn_param_traj->get_are_vals_fixed()) {
                        ixn_param_traj->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one_l2(timepoint,dt,options.l2_lambda,options.l2_center);
                    };
                };
            };
        } else {
            for (auto timepoint=timepoint_start_WSA + no_timesteps_WSA; timepoint>timepoint_start_WSA; timepoint--) {
                for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                    if (!ixn_param_traj->get_is_val_fixed_to_init_cond() && !ixn_param_traj->get_are_vals_fixed()) {
                        ixn_param_traj->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one(timepoint,dt);
                    };
                };
            };
        };
        
        if (options.verbose) {
            std::cout << "--- [Finished] ---" << std::endl;
        };

        /********************
         Form the update
         ********************/
        
        if (options.verbose) {
            std::cout << "--- Forming update ---" << std::endl;
        };

        for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
            if (!ixn_param_traj->get_is_val_fixed_to_init_cond() && !ixn_param_traj->get_are_vals_fixed()) {
                
                ixn_param_traj->get_diff_eq_rhs()->update_calculate_and_store(_timepoint_start_lattice,_timepoint_start_lattice+_no_timesteps_lattice,dt);
            };
        };
        
        if (options.verbose) {
            std::cout << "--- [Finished] ---" << std::endl;
        };

        /********************
         Committ the update
         ********************/
        
        if (options.verbose) {
            std::cout << "--- Committing update ---" << std::endl;
        };

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
        
        if (options.verbose) {
            std::cout << "--- [Finished] ---" << std::endl;
        };
    };
};
