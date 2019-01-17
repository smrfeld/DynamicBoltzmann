#include <string>
#include <vector>
#include <memory>
#include <map>

#ifndef SOLVER_H
#define SOLVER_H
#include "solver.hpp"
#endif

#ifndef OPT_PROBLEM_H
#define OPT_PROBLEM_H
#include "opt_problem.hpp"
#endif

/************************************
* Namespace for bmla
************************************/

namespace bmla {

	// Forwards
	class IxnParam;
	class Lattice;
    class FNameColl;

	/****************************************
	OptProblemRBM Options
	****************************************/

	struct OptionsWakeSleepRBM {

		// Verbosity
		bool verbose = false;
        
        // Is the hidden layer binary for the awake moment?
        bool is_awake_moment_hidden_binary = false;
        
        // Asleep
        CDModeAsleep cd_mode_asleep = CDModeAsleep::PERSISTENT_CD;
        
        // Sampling options
        // Is the visible reconstruction binary?
        bool is_asleep_visible_binary = true;
        // Is the hidden reconstruction binary, EXCEPT in the last phase?
        bool is_asleep_hidden_binary = true;
        // Is the hidden reconstruction binary in the last phase?
        bool is_asleep_hidden_binary_final = false;
        
        // Options for CD
        OptionsAsleepPersistentCD options_asleep_persistent_cd = OptionsAsleepPersistentCD();
        OptionsAsleepStartFromRandom options_asleep_start_from_random = OptionsAsleepStartFromRandom();
        OptionsAsleepStartFromData options_asleep_start_from_data = OptionsAsleepStartFromData();
        
        // Parallel
        bool parallel = true;
	};

	/****************************************
	OptProblemRBM
	****************************************/

    class OptProblemRBM : public OptProblem {

	public:

		/********************
		Constructor
		********************/

        using OptProblem::OptProblem;
        
		/********************
		Wake/asleep loop
		********************/

		void wake_sleep_loop(int batch_size, int no_markov_chains, int no_cd_sampling_steps, FNameColl &fname_coll, OptionsWakeSleepRBM options);

		/********************
		Solve
		********************/

		// Check if options passed are valid
		void check_options(int batch_size, int no_markov_chains, double dopt, int no_cd_sampling_steps, OptionsSolve options, OptionsWakeSleepRBM options_wake_sleep);

		// One step
		void solve_one_step(int i_opt_step, int batch_size, int no_markov_chains, double dopt, int no_cd_sampling_steps, FNameColl &fname_coll, OptionsSolve options, OptionsWakeSleepRBM options_wake_sleep);

		// Many steps
		void solve(int no_opt_steps, int batch_size, int no_markov_chains, double dopt, int no_cd_sampling_steps, FNameColl &fname_coll, OptionsSolve options, OptionsWakeSleepRBM options_wake_sleep);
	};

};
