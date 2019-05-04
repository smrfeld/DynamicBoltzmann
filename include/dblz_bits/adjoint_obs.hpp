#include <string>
#include <map>
#include <vector>

#include "fwds/fwds_ixn_param_traj.hpp"
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
    // MARK: - Common term
    // ***************
    
    class AdjointObsCommonTerm {
        
    private:
        
        // Species and index in the diff eq
        Sptr _species;
        int _idx_diff_eq;
        
        // All ixn param trajs
        std::vector<ITptr> _all_ixn_param_trajs;
        
        // Values
        // Timepoints = timesteps + 1
        int _no_timepoints;
        int _no_timesteps;
        std::vector<double> _vals;
        
        // Internal copy func/clean up
        void _clean_up();
        void _copy(const AdjointObsCommonTerm& other);
        void _move(AdjointObsCommonTerm &other);
        
    public:
        
        // ***************
        // MARK: - Constructor
        // ***************
        
        AdjointObsCommonTerm(Sptr species, int idx_diff_eq, std::vector<ITptr> all_ixn_param_trajs);
        AdjointObsCommonTerm(const AdjointObsCommonTerm& other);
        AdjointObsCommonTerm& operator=(const AdjointObsCommonTerm& other);
        AdjointObsCommonTerm(AdjointObsCommonTerm&& other);
        AdjointObsCommonTerm& operator=(AdjointObsCommonTerm&& other);
        ~AdjointObsCommonTerm();
        
        // ***************
        // MARK: - Timesteps
        // ***************
        
        int get_no_timesteps() const;
        void set_no_timesteps(int no_timesteps);
        
        // ***************
        // MARK: - Vals
        // ***************
        
        void calculate_val_at_timepoint(int timepoint);
        double get_val_at_timepoint(int timepoint) const;
        
    };
    
    // ***************
    // MARK: - Adjoint moment term
    // ***************
    
    class AdjointMomentTerm {
        
    private:
        
        // Ixn param traj - this is part of the sum
        ITptr _ixn_param;
        
        // The other part of the sum from the domain on the RHS of the diff eq
        int _layer;
        Sptr _species;
        
        // Vals
        int _no_timesteps;
        int _no_timepoints;
        std::vector<double> _vals;
        
        // Constructor helpers
        void _clean_up();
        void _move(AdjointMomentTerm &other);
        void _copy(const AdjointMomentTerm& other);
        
    public:
        
        /********************
         Constructor
         ********************/
        
        AdjointMomentTerm(ITptr ixn_param, int layer, Sptr species);
        AdjointMomentTerm(const AdjointMomentTerm& other);
        AdjointMomentTerm(AdjointMomentTerm&& other);
        AdjointMomentTerm& operator=(const AdjointMomentTerm& other);
        AdjointMomentTerm& operator=(AdjointMomentTerm&& other);
        ~AdjointMomentTerm();
        
        /********************
         Name, type
         ********************/
        
        ITptr get_ixn_param_traj() const;
        int get_layer() const;
        Sptr get_species() const;
        
        /********************
         Timepoints
         ********************/
        
        void set_no_timesteps(int no_timesteps);
        int get_no_timesteps() const;
        
        /********************
         Get/set moment
         ********************/
        
        // Get moment
        double get_val_at_timepoint(int timepoint) const;
        void increment_val_at_timepoint(int timepoint, double val);
        void set_val_at_timepoint(int timepoint, double val);
        void reset_val_at_timepoint(int timepoint);

    };
    
    // ***************
    // MARK: - Adjoint class when observables used for diff eq RHS
    // ***************
    
    class AdjointObs : public Adjoint {

	private:
        
        std::shared_ptr<AdjointObsCommonTerm> _common_term;

		// Internal copy func/clean up
		void _clean_up();
		void _copy(const AdjointObs& other);
		void _move(AdjointObs &other);

	public:

        // ***************
        // MARK: - Constructor
        // ***************

        AdjointObs(std::string name, ITptr ixn_param_traj, std::shared_ptr<AdjointObsCommonTerm> common_term);
		AdjointObs(const AdjointObs& other);
		AdjointObs& operator=(const AdjointObs& other);
		AdjointObs(AdjointObs&& other);
		AdjointObs& operator=(AdjointObs&& other);
		~AdjointObs();

        // ***************
        // MARK: - Solve diff eq
        // ***************

        void solve_diff_eq_at_timepoint_to_minus_one(int timepoint, double dt);
        void solve_diff_eq_at_timepoint_to_minus_one_l2(int timepoint, double dt, double l2_lambda, double l2_center);
	};
};
