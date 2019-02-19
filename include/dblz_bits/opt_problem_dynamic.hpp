#include <string>
#include <vector>
#include <memory>
#include <map>

/************************************
 * Namespace for bmla
 ************************************/

namespace dblz {
    
    // Forwards
    class IxnParamTraj;
    class LatticeTraj;
    class FNameTrajColl;
    struct OptionsWakeSleep;
    
    /****************************************
    Misc options
     ****************************************/

    enum class Solver : unsigned int { SGD, NESTEROV, ADAM };
    enum class MCType: unsigned int;
    
    /****************************************
    Options
     ****************************************/
    
    struct OptionsSolveDynamic {
        
        // Should check options before starting
        bool should_check_options = true;
        
        // Verbosity
        bool verbose = false;
        bool verbose_timing = true;
        
        // L2 Reg mode
        bool l2_reg = false;
        std::map<std::shared_ptr<IxnParamTraj>,double> l2_lambda;
        std::map<std::shared_ptr<IxnParamTraj>,double> l2_center;
        
        // Pointwise l2_reg
        bool l2_reg_center_traj = false;
        std::map<std::shared_ptr<IxnParamTraj>,std::map<int,double>> l2_center_traj;

        // Options for the solvers
        Solver solver = Solver::ADAM;
        
        // Nesterov
        double nesterov_acc = 0.5;
        
        // Adam
        double adam_beta_1 = 0.9;
        double adam_beta_2 = 0.999;
        double adam_eps = 0.00000001;
        
        // No steps per step for ixn params/adjoint
        int no_steps_per_step_IP = 1;
    };
    
    /****************************************
     OptProblemDynamic
     ****************************************/
    
    class OptProblemDynamic {
        
    protected:
        
        // Lattice
        std::shared_ptr<LatticeTraj> _latt_traj;
        
        // Constructor helpers
        void _clean_up();
        void _move(OptProblemDynamic &other);
        void _copy(const OptProblemDynamic& other);
        
    public:
        
        /********************
         Constructor
         ********************/
        
        OptProblemDynamic(std::shared_ptr<LatticeTraj> latt_traj);
        OptProblemDynamic(const OptProblemDynamic& other);
        OptProblemDynamic(OptProblemDynamic&& other);
        OptProblemDynamic& operator=(const OptProblemDynamic &other);
        OptProblemDynamic& operator=(OptProblemDynamic &&other);
        ~OptProblemDynamic();

        /********************
         Solve
         ********************/
        
        // Check if options passed are valid
        void check_options(int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_mean_field_updates, int no_gibbs_sampling_steps, OptionsSolveDynamic options, OptionsWakeSleep options_wake_sleep);
        
        // One step
        // SIP = solve ixn params
        // WSA = wake/sleep/adjoint
        void solve_one_step(int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_mean_field_updates, int no_gibbs_sampling_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep options_wake_sleep);
    };
    
};
