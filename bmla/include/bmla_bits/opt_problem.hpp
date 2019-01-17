#include <string>
#include <vector>
#include <memory>
#include <map>

#ifndef SOLVER_H
#define SOLVER_H
#include "solver.hpp"
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
     OptProblemRBM
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
    };
    
};
