#include <string>
#include <map>
#include <vector>

#include "fwds/fwds_ixn_param_traj.hpp"
#include "fwds/fwds_center_traj.hpp"
#include "fwds/fwds_species.hpp"

#ifndef ADJOINT_H
#define ADJOINT_H
#include "adjoint.hpp"
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {	

    // ***************
    // MARK: - Derivative term
    // ***************
    
    class AdjointParamsCenteredHomDerivTerm {
        
    private:
        
        // This ixn
        ITptr _deriv_ixn_param_traj;
        
        // All ixn params
        // Biases:
        // Layer -> species -> ixn
        std::map<int,std::map<Sptr,ITptr>> _all_biases;
        // Weights:
        // Layer l -> species layer l -> species layer l+1 -> ixn and idx of deriv
        // weight l->l+1
        std::map<int, std::map<Sptr, std::map<int, std::map<Sptr,ITptr>>>> _all_weights;

        // All centers
        // Layer -> species -> center
        std::map<int,std::map<Sptr,CTptr>> _all_center_trajs;
        
        // Values
        // Timepoints = timesteps + 1
        int _no_timepoints;
        int _no_timesteps;
        std::vector<double> _vals;
        
        // Connection multiplicities across layers
        // layer l -> layer l+1 or l-1 -> mult
        std::map<int, std::map<int,int>> _conn_mults;
        
        // Internal copy func/clean up
        void _clean_up();
        void _copy(const AdjointParamsCenteredHomDerivTerm& other);
        void _move(AdjointParamsCenteredHomDerivTerm &other);
        
    public:
        
        // ***************
        // MARK: - Constructor
        // ***************
        
        AdjointParamsCenteredHomDerivTerm(ITptr deriv_ixn_param_traj, std::map<int,std::map<Sptr,ITptr>> all_biases, std::map<int, std::map<Sptr, std::map<int, std::map<Sptr,ITptr>>>> all_weights, std::map<int,std::map<Sptr,CTptr>> all_center_trajs, std::map<int, std::map<int,int>> conn_mults);
        AdjointParamsCenteredHomDerivTerm(const AdjointParamsCenteredHomDerivTerm& other);
        AdjointParamsCenteredHomDerivTerm& operator=(const AdjointParamsCenteredHomDerivTerm& other);
        AdjointParamsCenteredHomDerivTerm(AdjointParamsCenteredHomDerivTerm&& other);
        AdjointParamsCenteredHomDerivTerm& operator=(AdjointParamsCenteredHomDerivTerm&& other);
        AdjointParamsCenteredHomDerivTerm();
        
        // ***************
        // MARK: - Getters
        // ***************
        
        ITptr get_deriv_ixn_param_traj() const;
        
        // ***************
        // MARK: - Timesteps
        // ***************
        
        int get_no_timesteps() const;
        void set_no_timesteps(int no_timesteps);
        
        // ***************
        // MARK: - Vals
        // ***************
        
        void calculate_val_at_timepoint(int timepoint, bool form_abscissas=true);
        double get_val_at_timepoint(int timepoint) const;
    };
    
    // ***************
    // MARK: - Adjoint class when params used for diff eq RHS
    // ***************
    
    typedef std::shared_ptr<AdjointParamsCenteredHomDerivTerm> CtrDerivPtr;
    
    class AdjointParamsCenteredHomBias : public Adjoint {

	private:
        
        CtrDerivPtr _deriv_term_bias;

		// Internal copy func/clean up
		void _clean_up();
		void _copy(const AdjointParamsCenteredHomBias& other);
		void _move(AdjointParamsCenteredHomBias &other);

	public:

        // ***************
        // MARK: - Constructor
        // ***************

		AdjointParamsCenteredHomBias(std::string name, ITptr ixn_param_traj, CtrDerivPtr deriv_term_bias);
		AdjointParamsCenteredHomBias(const AdjointParamsCenteredHomBias& other);
		AdjointParamsCenteredHomBias& operator=(const AdjointParamsCenteredHomBias& other);
		AdjointParamsCenteredHomBias(AdjointParamsCenteredHomBias&& other);
		AdjointParamsCenteredHomBias& operator=(AdjointParamsCenteredHomBias&& other);
		~AdjointParamsCenteredHomBias();

        // ***************
        // MARK: - Get deriv term
        // ***************
        
        CtrDerivPtr get_deriv_term_bias() const;
        
        // ***************
        // MARK: - Solve diff eq
        // ***************

		void solve_diff_eq_at_timepoint_to_minus_one(int timepoint, double dt, bool form_abscissas=true);
        void solve_diff_eq_at_timepoint_to_minus_one_l2(int timepoint, double dt, double l2_lambda, double l2_center, bool form_abscissas=true);
	};
    
    class AdjointParamsCenteredHomWeight : public Adjoint {
        
    private:
        
        CtrDerivPtr _deriv_term_weight;
        CtrDerivPtr _deriv_term_bias_lower;
        CtrDerivPtr _deriv_term_bias_upper;
        int _conn_mult;
        CTptr _center_lower;
        CTptr _center_upper;
        std::shared_ptr<AdjointParamsCenteredHomBias> _adjoint_bias_lower;
        std::shared_ptr<AdjointParamsCenteredHomBias> _adjoint_bias_upper;
        
        // Internal copy func/clean up
        void _clean_up();
        void _copy(const AdjointParamsCenteredHomWeight& other);
        void _move(AdjointParamsCenteredHomWeight &other);
        
    public:
        
        // ***************
        // MARK: - Constructor
        // ***************
        
        AdjointParamsCenteredHomWeight(std::string name, ITptr ixn_param_traj, CtrDerivPtr deriv_term_weight, CtrDerivPtr deriv_term_bias_lower, CtrDerivPtr deriv_term_bias_upper, int conn_mult, CTptr center_lower, CTptr center_upper, std::shared_ptr<AdjointParamsCenteredHomBias> adjoint_bias_lower, std::shared_ptr<AdjointParamsCenteredHomBias> adjoint_bias_upper);
        AdjointParamsCenteredHomWeight(const AdjointParamsCenteredHomWeight& other);
        AdjointParamsCenteredHomWeight& operator=(const AdjointParamsCenteredHomWeight& other);
        AdjointParamsCenteredHomWeight(AdjointParamsCenteredHomWeight&& other);
        AdjointParamsCenteredHomWeight& operator=(AdjointParamsCenteredHomWeight&& other);
        ~AdjointParamsCenteredHomWeight();
        
        // ***************
        // MARK: - Get deriv term
        // ***************
        
        CtrDerivPtr get_deriv_term_weight() const;
        CtrDerivPtr get_deriv_term_bias_lower() const;
        CtrDerivPtr get_deriv_term_bias_upper() const;

        // ***************
        // MARK: - Solve diff eq
        // ***************
        
        void solve_diff_eq_at_timepoint_to_minus_one(int timepoint, double dt, bool form_abscissas=true);
        void solve_diff_eq_at_timepoint_to_minus_one_l2(int timepoint, double dt, double l2_lambda, double l2_center, bool form_abscissas=true);
    };
};
