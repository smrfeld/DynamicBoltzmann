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
    class LatticeTraj1DFullyVisible;
    class FNameTrajColl;
    struct OptionsWakeSleep_BM_PCD;
    struct OptionsWakeSleep_RBM_CD;
    struct OptionsWakeSleep_1DFV_CD;

    /****************************************
    Misc options
     ****************************************/

    enum class Solver : unsigned int { SGD, ADAM };
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
        
        // L2 traj
        bool l2_reg_traj = false;
        std::map<std::shared_ptr<IxnParamTraj>,std::map<int,double>> l2_lambda_traj;
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
        
        // Locking mode is on/off
        bool locking_mode = false;
        
        // Sliding factor
        double sliding_factor = 0.01;
    };
    
    /****************************************
     OptProblemDynamic
     ****************************************/
    
    class OptProblemDynamic {
        
    protected:
        
        // Lattice
        std::shared_ptr<LatticeTraj> _latt_traj;
        std::shared_ptr<LatticeTraj1DFullyVisible> _latt_traj_1dfv;

        // Constructor helpers
        void _clean_up();
        void _move(OptProblemDynamic &other);
        void _copy(const OptProblemDynamic& other);
        
    public:
        
        // ***************
        // MARK: - Constructor
        // ***************
      
        OptProblemDynamic(std::shared_ptr<LatticeTraj> latt_traj);
        OptProblemDynamic(std::shared_ptr<LatticeTraj1DFullyVisible> latt_traj_1dfv);
        OptProblemDynamic(const OptProblemDynamic& other);
        OptProblemDynamic(OptProblemDynamic&& other);
        OptProblemDynamic& operator=(const OptProblemDynamic &other);
        OptProblemDynamic& operator=(OptProblemDynamic &&other);
        ~OptProblemDynamic();

        // ***************
        // MARK: - Solve ixn params
        // ***************
        
        // Solve t to t+delta t
        void solve_ixn_param_trajs_step(const std::vector<std::shared_ptr<IxnParamTraj>> &ixn_param_trajs, double dt, int timepoint, int no_steps_per_step) const;
        // Solve many points
        void solve_ixn_param_trajs(const std::vector<std::shared_ptr<IxnParamTraj>> &ixn_param_trajs, double dt, int timepoint_start, int no_timesteps, int no_steps_per_step) const;

        // ***************
        // MARK: - Committ step
        // ***************
        
        void committ_step(const std::vector<std::shared_ptr<IxnParamTraj>> &ixn_params, int i_opt_step, OptionsSolveDynamic options);

        // ***************
        // MARK: - BM PCD params
        // ***************
        
        // One step
        // SIP = solve ixn params
        // WSA = wake/sleep/adjoint
        void solve_one_step_bm_pcd_params(int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_mean_field_updates, int no_gibbs_sampling_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_BM_PCD options_wake_sleep);
        void solve_one_step_bm_pcd_params_without_committ(int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_mean_field_updates, int no_gibbs_sampling_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_BM_PCD options_wake_sleep);
        
        // ***************
        // MARK: - RBM CD params
        // ***************
        
        void solve_one_step_rbm_cd_params(int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_cd_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_RBM_CD options_wake_sleep);
        void solve_one_step_rbm_cd_params_without_committ(int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_cd_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_RBM_CD options_wake_sleep);
        
        // ***************
        // MARK: - RBM CD obs
        // ***************
        
        void solve_one_step_rbm_cd_obs(int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_cd_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_RBM_CD options_wake_sleep);
        void solve_one_step_rbm_cd_obs_without_committ(int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_cd_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_RBM_CD options_wake_sleep);
        
        // ***************
        // MARK: - 1D fully visible lattice
        // ***************
        
        void solve_one_step_1d_fully_visible(int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_cd_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_1DFV_CD options_wake_sleep);
        void solve_one_step_1d_fully_visible_without_committ(int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_cd_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_1DFV_CD options_wake_sleep);
        
        // ***************
        // MARK: - RBM CD Centered Hom
        // ***************
        
        void solve_one_step_rbm_cd_centered_hom(int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_cd_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_RBM_CD options_wake_sleep);
        void solve_one_step_rbm_cd_centered_hom_without_committ(int i_opt_step, int timepoint_start_SIP, int no_timesteps_SIP, int timepoint_start_WS, int no_timesteps_WS, int timepoint_start_A, int no_timesteps_A, double dt, int no_cd_steps, FNameTrajColl &fname_traj_coll, OptionsSolveDynamic options, OptionsWakeSleep_RBM_CD options_wake_sleep);
    };
    
};
