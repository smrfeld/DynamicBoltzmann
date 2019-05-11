#include <vector>
#include <string>

#include "fwds/fwds_species.hpp"
#include "fwds/fwds_ixn_param.hpp"

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	// Forwards
	class DiffEqRHS;
	class Adjoint;
    class AdjointObs;
    class AdjointParams;
    class MomentDiff;
    
	/****************************************
	Interaction parameter
	****************************************/

	// Enumeration of type of dimension
    enum class IxnParamType: unsigned int;

	class IxnParamTraj {

	private:
        
        // Adjoint
        std::shared_ptr<Adjoint> _adjoint;
        std::shared_ptr<AdjointObs> _adjoint_obs;
        std::shared_ptr<AdjointParams> _adjoint_params;

        // Diff eq RHS
        std::shared_ptr<DiffEqRHS> _diff_eq;
        
        // Diff eq rhs that this ixn param appears in
        // and the corresponding idx in that function (= 0, 1, 2, etc to denote argument order)
        std::vector<std::pair<std::shared_ptr<DiffEqRHS>,int>> _diff_eq_dependencies;

        // Ixn params
        std::vector<std::shared_ptr<IxnParam>> _ixn_params;
        
        // Intermediate substep vals
        double _substep_val_current;
        double _substep_val_new;
        
        // Init cond
        double _init_cond;
        
        // No timesteps/timepoints
        int _no_timesteps;
        int _no_timepoints; // = timepoints + 1
        
        // Fixed value
        bool _is_val_fixed;

        // Copy, clean up
        void _clean_up();
        void _copy(const IxnParamTraj& other);
        void _move(IxnParamTraj &other);

	public:

        // ***************
        // MARK: - Constructor
        // ***************
        
		IxnParamTraj(std::string name, IxnParamType type, double init_cond);
		IxnParamTraj(const IxnParamTraj& other);
		IxnParamTraj(IxnParamTraj&& other);
		IxnParamTraj& operator=(const IxnParamTraj& other);
		IxnParamTraj& operator=(IxnParamTraj&& other);
		~IxnParamTraj();
		
        // ***************
        // MARK: - Print moment string
        // ***************
        
        void print_val_traj(int timepoint_start, int no_timesteps) const;
        void print_moment_traj(int timepoint_start, int no_timesteps) const;
        void print_moment_diff_traj(int timepoint_start, int no_timesteps) const;

        // ***************
        // MARK: - Substep vals
        // ***************
        
        double get_substep_val() const;
        void set_substep_val_new(double val);
        void set_substep_val_current(double val);
        void move_to_new_substep_val();
        
        // ***************
        // MARK: - Diff eq rhs that this ixn param appears in
        // ***************

		void add_diff_eq_dependency(std::shared_ptr<DiffEqRHS> diff_eq, int idx);
		const std::vector<std::pair<std::shared_ptr<DiffEqRHS>,int>>& get_diff_eq_dependencies() const;

        // ***************
        // MARK: - Fix
        // ***************
        
        void set_fix_value(bool fixed);
        bool get_is_val_fixed() const;
        
        // ***************
        // MARK: - Timesteps
        // ***************
        
		int get_no_timesteps() const;
		void set_no_timesteps(int no_timesteps);

        // ***************
        // MARK: - Init cond
        // ***************

		double get_init_cond() const;
		void set_init_cond(double init_cond);

        // ***************
        // MARK: - Name, type
        // ***************

		std::string get_name() const;

		IxnParamType get_type() const;

        // ***************
        // MARK: - Diff eq
        // ***************
        
		std::shared_ptr<DiffEqRHS> get_diff_eq_rhs() const;
		void set_diff_eq_rhs(std::shared_ptr<DiffEqRHS> diff_eq);

		void solve_diff_eq_at_timepoint_to_plus_one(int timepoint, double dt);

        // ***************
        // MARK: - Ixn param
        // ***************

        std::shared_ptr<IxnParam> get_ixn_param_at_timepoint(int timepoint) const;
        
        // ***************
        // MARK: - Adjoint
        // ***************
        
		void set_adjoint(std::shared_ptr<AdjointObs> adjoint_obs);
        void set_adjoint(std::shared_ptr<AdjointParams> adjoint_params);
        std::shared_ptr<Adjoint> get_adjoint() const;
        std::shared_ptr<AdjointObs> get_adjoint_obs() const;
        std::shared_ptr<AdjointParams> get_adjoint_params() const;

        // ***************
        // MARK: - Write to file
        // ***************

		void write_val_traj_to_file(int timepoint_start, int no_timesteps, std::string fname, bool with_timepoints=true) const;
        void write_moment_traj_to_file(int timepoint_start, int no_timesteps, std::string fname, bool with_timepoints=true) const;
        void write_adjoint_traj_to_file(int timepoint_start, int no_timesteps, std::string fname, bool with_timepoints=true) const;
	};

};


