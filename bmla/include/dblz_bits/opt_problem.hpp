#include <string>
#include <vector>
#include <memory>
#include <map>

/************************************
 * Namespace for bmla
 ************************************/

namespace dblz {
    
    // Forwards
    class IxnParam;
    class Lattice;
    class FNameColl;
    
    /****************************************
    Misc options
     ****************************************/

    enum class Solver : unsigned int { SGD, NESTEROV, ADAM };
    enum class MCType: unsigned int;
    
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
        
        // Sampling options
        // Is the visible reconstruction binary, EXCEPT in the last phase?
        bool is_asleep_visible_binary = true;
        // Is the hidden reconstruction binary, EXCEPT in the last phase?
        bool is_asleep_hidden_binary = true;
        // Is the visible reconstruction binary in the last phase?
        bool is_asleep_visible_binary_final = true;
        // Is the hidden reconstruction binary in the last phase?
        bool is_asleep_hidden_binary_final = false;

        // Write after awake/asleep
        bool write_after_awake = false;
        bool write_after_asleep = false;
        std::string write_after_awake_dir = "";
        std::string write_after_asleep_dir = "";
    };

    /****************************************
     OptProblem
     ****************************************/
    
    class OptProblem {
        
    protected:
        
        // No markov chains
        std::map<MCType,int> _no_markov_chains;
        
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
        
        OptProblem(std::shared_ptr<Lattice> latt, int no_markov_chains_awake, int no_markov_chains_asleep);
        OptProblem(const OptProblem& other);
        OptProblem(OptProblem&& other);
        OptProblem& operator=(const OptProblem &other);
        OptProblem& operator=(OptProblem &&other);
        ~OptProblem();
        
        /********************
         Init structures
         ********************/
        
        void set_no_markov_chains(MCType chain, int no_markov_chains);
        
        /********************
         Wake/asleep loop
         ********************/
        
        void wake_sleep_loop(int i_opt_step, int no_mean_field_updates, int no_gibbs_sampling_steps, FNameColl &fname_coll, OptionsWakeSleep options);
        
        /********************
         Solve
         ********************/
        
        // Check if options passed are valid
        void check_options(double dopt, int no_mean_field_updates, int no_gibbs_sampling_steps, OptionsSolve options, OptionsWakeSleep options_wake_sleep);
        
        // One step
        void solve_one_step(int i_opt_step, double dopt, int no_mean_field_updates, int no_gibbs_sampling_steps, FNameColl &fname_coll, OptionsSolve options, OptionsWakeSleep options_wake_sleep);
        
        // Many steps
        void solve(int no_opt_steps, double dopt, int no_mean_field_updates, int no_gibbs_sampling_steps, FNameColl &fname_coll, OptionsSolve options, OptionsWakeSleep options_wake_sleep);
    };
    
};
