#include <string>
#include <vector>
#include <memory>
#include <map>

/************************************
 * Namespace for bmla
 ************************************/

namespace bmla {
    
    // Forwards
    class IxnParam;
    class Lattice;
    class FNameColl;
    
    /****************************************
    Misc options
     ****************************************/

    enum class Solver : unsigned int { SGD, NESTEROV, ADAM };
    
    enum class CDModeAsleep : unsigned int { PERSISTENT_CD, START_FROM_DATA, START_FROM_RANDOM };
    
    struct OptionsAsleepPersistentCD {};
    
    struct OptionsAsleepStartFromRandom {
        
        // Init random - binary?
        bool start_from_binary_visible = true;
        // Binary hidden?
        bool start_from_binary_hidden = true;
    };
    
    struct OptionsAsleepStartFromData {
        // Start from binary?
        bool start_from_binary_hidden = true;
    };
    
    /****************************************
    Options
     ****************************************/
    
    struct OptionsSolve {
        
        // Should check options before starting
        bool should_check_options = true;
        
        // Verbosity
        bool verbose_update = false;
        bool verbose_moment = true;
        
        // L2 Reg mode
        bool l2_reg = false;
        std::map<std::shared_ptr<IxnParam>,double> l2_lambda;
        std::map<std::shared_ptr<IxnParam>,double> l2_center;
        
        // Variable learning rate
        bool var_learning_rates = false;
        std::map<std::shared_ptr<IxnParam>,double> var_learning_rate_values;
        
        // Options for the solvers
        Solver solver = Solver::ADAM;
        
        // Nesterov
        double nesterov_acc = 0.5;
        
        // Adam
        double adam_beta_1 = 0.9;
        double adam_beta_2 = 0.999;
        double adam_eps = 0.00000001;
    };
    
    /****************************************
     OptProblem Options
     ****************************************/
    
    struct OptionsWakeSleep {
        
        // Verbosity
        bool verbose = false;
        
        // Mode
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
    };

    /****************************************
     OptProblem
     ****************************************/
    
    class OptProblem {
        
    protected:
        
        // Ixn params
        std::vector<std::shared_ptr<IxnParam>> _ixn_params;
        
        // Lattice
        std::shared_ptr<Lattice> _latt;
        
        // Constructor helpers
        void _clean_up();
        void _move(OptProblem &other);
        void _copy(const OptProblem& other);
        
    public:
        
        /********************
         Constructor
         ********************/
        
        OptProblem(std::shared_ptr<Lattice> latt, std::vector<std::shared_ptr<IxnParam>> ixn_params);
        OptProblem(const OptProblem& other);
        OptProblem(OptProblem&& other);
        OptProblem& operator=(const OptProblem &other);
        OptProblem& operator=(OptProblem &&other);
        ~OptProblem();
        
        /********************
         Init structures
         ********************/
        
        void init_structures(int batch_size, int no_markov_chains);
        
        /********************
         Wake/asleep loop
         ********************/
        
        void wake_sleep_loop(int batch_size, int no_markov_chains, int no_cd_sampling_steps, int no_mean_field_updates, FNameColl &fname_coll, OptionsWakeSleep options);
        
        /********************
         Solve
         ********************/
        
        // Check if options passed are valid
        void check_options(int batch_size, int no_markov_chains, double dopt, int no_cd_sampling_steps, int no_mean_field_updates, OptionsSolve options, OptionsWakeSleep options_wake_sleep);
        
        // One step
        void solve_one_step(int i_opt_step, int batch_size, int no_markov_chains, double dopt, int no_cd_sampling_steps, int no_mean_field_updates, FNameColl &fname_coll, OptionsSolve options, OptionsWakeSleep options_wake_sleep);
        
        // Many steps
        void solve(int no_opt_steps, int batch_size, int no_markov_chains, double dopt, int no_cd_sampling_steps, int no_mean_field_updates, FNameColl &fname_coll, OptionsSolve options, OptionsWakeSleep options_wake_sleep);
    };
    
};
