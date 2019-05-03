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
    class Lattice1DFullyVisible;
    class LatticeCenteredHom;
    class FNameColl;
    struct OptionsWakeSleep_BM_PCD;
    struct OptionsWakeSleep_RBM_CD;
    struct OptionsWakeSleep_1DFV_CD;
    struct OptionsWakeSleep_RBM_CD_CH;

    /****************************************
    Misc options
     ****************************************/

    enum class Solver : unsigned int { SGD, ADAM };
    enum class MCType: unsigned int;
    
    /****************************************
    Options
     ****************************************/
    
    struct OptionsSolveStatic {
        
        // Should check options before starting
        bool should_check_options = true;
        
        // Verbosity
        bool verbose_update = false;
        bool verbose_moment = true;
        
        // L2 Reg mode
        bool l2_reg = false;
        std::map<std::shared_ptr<IxnParam>,double> l2_lambda;
        std::map<std::shared_ptr<IxnParam>,double> l2_center;
                
        // Options for the solvers
        Solver solver = Solver::ADAM;
        
        // Nesterov
        double nesterov_acc = 0.5;
        
        // Adam
        double adam_beta_1 = 0.9;
        double adam_beta_2 = 0.999;
        double adam_eps = 0.00000001;
        
        // Sliding factor
        double sliding_factor = 0.01;
    };
    
    /****************************************
     OptProblemStatic
     ****************************************/
    
    class OptProblemStatic {
        
    protected:
        
        // Lattice
        std::shared_ptr<Lattice> _latt;
        std::shared_ptr<Lattice1DFullyVisible> _latt1dfv;
        std::shared_ptr<LatticeCenteredHom> _lattch;

        // Constructor helpers
        void _clean_up();
        void _move(OptProblemStatic &other);
        void _copy(const OptProblemStatic& other);
        
    public:
        
        /********************
         Constructor
         ********************/
        
        OptProblemStatic(std::shared_ptr<Lattice> latt);
        OptProblemStatic(std::shared_ptr<Lattice1DFullyVisible> latt1dfv);
        OptProblemStatic(std::shared_ptr<LatticeCenteredHom> lattch);
        OptProblemStatic(const OptProblemStatic& other);
        OptProblemStatic(OptProblemStatic&& other);
        OptProblemStatic& operator=(const OptProblemStatic &other);
        OptProblemStatic& operator=(OptProblemStatic &&other);
        ~OptProblemStatic();
        
        /********************
         Solve
         ********************/
        
        // One step
        void solve_one_step_bm_pcd(int i_opt_step, int no_mean_field_updates, int no_gibbs_sampling_steps, FNameColl &fname_coll, OptionsSolveStatic options, OptionsWakeSleep_BM_PCD options_wake_sleep);

        void solve_one_step_rbm_cd(int i_opt_step, int no_cd_steps, FNameColl &fname_coll, OptionsSolveStatic options, OptionsWakeSleep_RBM_CD options_wake_sleep);
        
        void solve_one_step_1d_fully_visible(int i_opt_step, int no_cd_steps, FNameColl &fname_coll, OptionsSolveStatic options, OptionsWakeSleep_1DFV_CD options_wake_sleep);

        void solve_one_step_rbm_cd_centered_hom(int i_opt_step, int no_cd_steps, FNameColl &fname_coll, OptionsSolveStatic options, OptionsWakeSleep_RBM_CD_CH options_wake_sleep);
    };
};
