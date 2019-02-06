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
        
        set_no_markov_chains(MCType::AWAKE,no_markov_chains_awake);
        set_no_markov_chains(MCType::ASLEEP,no_markov_chains_asleep);
        set_no_timesteps_ixn_params(no_timesteps_ixn_params);
        set_no_timesteps_lattice(timepoint_start_lattice, no_timesteps_lattice);
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
        for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
            for (auto timepoint=0; timepoint<_no_timepoints_ixn_params; timepoint++) {
                ixn_param_traj->get_ixn_param_at_timepoint(timepoint)->get_moment()->set_no_markov_chains(chain, no_markov_chains);
            };
        };
        
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
     Wake/asleep loop
     ********************/
    
    void OptProblemDynamic::wake_sleep_loop(int i_opt_step, int no_mean_field_updates, int no_gibbs_sampling_steps, FNameTrajColl &fname_traj_coll, OptionsWakeSleepDynamic options) {
        
        // Reset all moments
        for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
            for (auto timepoint=_timepoint_start_lattice; timepoint<_timepoint_start_lattice+_no_timesteps_lattice; timepoint++) {
                auto moment = ixn_param_traj->get_ixn_param_at_timepoint(timepoint)->get_moment();
                moment->reset_to_zero(MCType::ASLEEP);
                if (!moment->get_is_awake_moment_fixed()) {
                    moment->reset_to_zero(MCType::AWAKE);
                };
            };
        };
        
        // AWAKE PHASE
        
        // Make a batch subset
        std::vector<int> idx_subset = fname_traj_coll.get_random_subset(_no_markov_chains[MCType::AWAKE]);
        
        // Go through all timepoints
        for (auto timepoint=_timepoint_start_lattice; timepoint<_timepoint_start_lattice+_no_timesteps_lattice; timepoint++) {

            auto latt = _latt_traj->get_lattice_at_timepoint(timepoint);
            
            // Read in the batch at all timepoints
            for (int i_chain=0; i_chain<_no_markov_chains[MCType::AWAKE]; i_chain++)
            {
                auto file_traj = fname_traj_coll.get_fname_traj(idx_subset[i_chain]);
                latt->read_layer_from_file(MCType::AWAKE, i_chain, 0, file_traj[timepoint].name, file_traj[timepoint].binary);
            };
            
            // Option (1): init MF with random hidden layers with prob units
            /*
            for (int i_chain=0; i_chain<_no_markov_chains[MCType::AWAKE]; i_chain++) {
                _latt_traj->set_random_all_hidden_units(MCType::AWAKE, i_chain, false);
            };
             */
            // Option (2): upward pass with 2x weights (DBM) to activate probabilitsic units
            // (faster to converge!!!)
            latt->activate_upward_pass_with_2x_weights_1x_bias(MCType::AWAKE, false);
            
            // Variational inference
            for (auto i=0; i<no_mean_field_updates; i++) {
                latt->mean_field_hiddens_step();
            };
            
            // Reap awake
            latt->reap_moments(MCType::AWAKE);
            
            // ASLEEP PHASE - PERSISTENT_CD
            
            // Run CD sampling
            for (auto i_chain=0; i_chain<_no_markov_chains[MCType::ASLEEP]; i_chain++) {
                
                // Sample vis, hidden
                for (int i_sampling_step=0; i_sampling_step<no_gibbs_sampling_steps; i_sampling_step++)
                {
                    if (i_sampling_step != no_gibbs_sampling_steps-1) {
                        latt->gibbs_sampling_step(options.is_asleep_visible_binary, options.is_asleep_hidden_binary);
                    } else {
                        if (options.is_asleep_visible_binary_final && options.is_asleep_hidden_binary_final) {
                            // All binary
                            latt->gibbs_sampling_step(options.is_asleep_visible_binary_final, options.is_asleep_hidden_binary_final);
                        } else {
                            // Parallel for non-binary options
                            latt->gibbs_sampling_step_parallel(options.is_asleep_visible_binary_final, options.is_asleep_hidden_binary_final);
                        };
                    };
                };
            };
            
            // Reap asleep
            latt->reap_moments(MCType::ASLEEP);
        
            // Average moments
            for (auto &ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                // std::cout << "averaging: " << ixn_param->get_name() << std::endl;
                if (!ixn_param_traj->get_ixn_param_at_timepoint(timepoint)->get_moment()->get_is_awake_moment_fixed()) {
                    ixn_param_traj->get_ixn_param_at_timepoint(timepoint)->get_moment()->average_moment_samples(MCType::AWAKE);
                };
                ixn_param_traj->get_ixn_param_at_timepoint(timepoint)->get_moment()->average_moment_samples(MCType::ASLEEP);
            };
        };
    };
    
    /********************
     Solve
     ********************/
    
    // Check if options passed are valid
    void OptProblemDynamic::check_options(double dopt, double dt, int no_mean_field_updates, int no_gibbs_sampling_steps, OptionsSolveDynamic options, OptionsWakeSleepDynamic options_wake_sleep) {
        if (options.solver == Solver::SGD || options.solver == Solver::NESTEROV) {
            std::cerr << ">>> OptProblemDynamic::check_options <<< Error: only Adam is currently supported as a solver" << std::endl;
            exit(EXIT_FAILURE);
        };
    };
    
    // One step
    void OptProblemDynamic::solve_one_step(int i_opt_step, double dt, double dopt, int no_mean_field_updates, int no_gibbs_sampling_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleepDynamic options_wake_sleep) {
        
        /*****
         Check options
         *****/
        
        if (options.should_check_options) {
            check_options(dopt,dt,no_mean_field_updates,no_gibbs_sampling_steps,options,options_wake_sleep);
        };
        
        /*****
         Solve diff eq for F
         *****/
        
        // Solve over all time
        for (auto timepoint=0; timepoint<_no_timesteps_ixn_params; timepoint++) {
            for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                if (!ixn_param_traj->get_is_val_fixed_to_init_cond() && !ixn_param_traj->get_are_vals_fixed()) {
                    ixn_param_traj->solve_diff_eq_at_timepoint_to_plus_one(timepoint,dt);
                };
            };
        };
        
        /*****
         Wake/asleep loop
         *****/
        
        wake_sleep_loop(i_opt_step,no_mean_field_updates,no_gibbs_sampling_steps,fname_traj_coll,options_wake_sleep);
        
        /********************
         Solve diff eq for adjoint
         ********************/
        
        if (options.l2_reg) {
            for (auto timepoint=_timepoint_start_lattice+_no_timesteps_lattice; timepoint>_timepoint_start_lattice; timepoint--) {
                for (auto ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                    if (!ixn_param_traj->get_is_val_fixed_to_init_cond() && !ixn_param_traj->get_are_vals_fixed()) {
                        ixn_param_traj->get_adjoint()->solve_diff_eq_at_timepoint_to_minus_one(timepoint,dt,options.l2_reg,options.l2_lambda,options.l2_center);
                    };
                };
            };
        } else {
            for (auto timepoint=_timepoint_start_lattice+_no_timesteps_lattice; timepoint>_timepoint_start_lattice; timepoint--) {
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
                
                ixn_param_traj->get_diff_eq_rhs()->update_calculate_and_store(_timepoint_start_lattice,_timepoint_start_lattice+_no_timepoints_lattice,dt);
            };
        };

        /********************
         Committ the update
         ********************/
        
        double dopt_use = dopt;
        if (options.solver == Solver::ADAM) {
            for (auto &ixn_param_traj: _latt_traj->get_all_ixn_param_trajs()) {
                if (!ixn_param_traj->get_is_val_fixed_to_init_cond() && !ixn_param_traj->get_are_vals_fixed()) {
                    // Learning rate
                    if (options.var_learning_rates) {
                        dopt_use = options.var_learning_rate_values[ixn_param_traj];
                    };
                    
                    ixn_param_traj->get_diff_eq_rhs()->update_committ_stored_adam(dopt_use,i_opt_step,options.adam_beta_1,options.adam_beta_2,options.adam_eps);
                };
            };
        };
    };
    
    // Many steps
    void OptProblemDynamic::solve(int no_opt_steps, double dt, double dopt, int no_mean_field_updates, int no_gibbs_sampling_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleepDynamic options_wake_sleep) {
        
        for (int i_opt_step=1; i_opt_step<=no_opt_steps; i_opt_step++)
        {
            
            std::cout << "------------------" << std::endl;
            std::cout << "Opt step: " << i_opt_step << " / " << no_opt_steps << std::endl;
            std::cout << "------------------" << std::endl;
            
            // Solve
            solve_one_step(i_opt_step,dt,dopt,no_mean_field_updates,no_gibbs_sampling_steps,fname_traj_coll,options,options_wake_sleep);
        };
    };
};
