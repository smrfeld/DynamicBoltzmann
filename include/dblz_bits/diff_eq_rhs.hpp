#include <vector>
#include <string>
#include <map>

#include "fwds/fwds_ixn_param_traj.hpp"
#include "fwds/fwds_species.hpp"
#include "fwds/fwds_center_traj.hpp"

#include <q3c1> // guard is built-in

/************************************
* Namespace for dblz
************************************/

namespace dblz {

    /****************************************
     Domain1D - a abstract base
     ****************************************/
    
    class Domain1D : public q3c1::Dimension1D {
        
    public:
        
        /********************
         Constructor
         ********************/
        
        // Inherit!
        Domain1D(double delta, double zero);
        virtual ~Domain1D();
        
        /********************
         Getters
         ********************/
        
        virtual double get_val_at_timepoint(int timepoint) const = 0;
    };

    
	/****************************************
	Domain1DParam
	****************************************/

	class Domain1DParam : public Domain1D {

	private:

		ITptr _ixn_param_traj;

		// Internal copy func/clean up
		void _clean_up();
		void _move(Domain1DParam& other);
		void _copy(const Domain1DParam& other);

	public:

		/********************
		Constructor
		********************/

		Domain1DParam(ITptr ixn_param_traj, double delta, double zero);
		Domain1DParam(const Domain1DParam& other);
		Domain1DParam& operator=(const Domain1DParam& other);
		Domain1DParam(Domain1DParam&& other);
		Domain1DParam& operator=(Domain1DParam&& other);
		~Domain1DParam();

		/********************
		Getters
		********************/

		ITptr get_ixn_param_traj() const;
        
        double get_val_at_timepoint(int timepoint) const;
	};

    /****************************************
     Domain1DObs
     ****************************************/
    
    class Domain1DObs : public Domain1D {
        
    private:
        
        // For now, only "sum of single species" moment is allowed
        int _layer;
        Sptr _species;
        
        // Values
        int _no_timesteps;
        int _no_timepoints;
        std::vector<double> _vals;
        
        // Internal copy func/clean up
        void _clean_up();
        void _move(Domain1DObs& other);
        void _copy(const Domain1DObs& other);
        
    public:
        
        /********************
         Constructor
         ********************/
        
        Domain1DObs(int layer, Sptr species, double val_at_time_zero, double delta, double zero);
        Domain1DObs(const Domain1DObs& other);
        Domain1DObs& operator=(const Domain1DObs& other);
        Domain1DObs(Domain1DObs&& other);
        Domain1DObs& operator=(Domain1DObs&& other);
        ~Domain1DObs();
        
        /********************
         Timesteps
         ********************/
        
        void set_no_timesteps(int no_timesteps);
        int get_no_timesteps() const;
        
        /********************
         Getters
         ********************/
        
        int get_layer() const;
        Sptr get_species() const;
        
        void set_val_at_timepoint(int timepoint, double val);
        double get_val_at_timepoint(int timepoint) const;
    };

    /****************************************
     Domain1DCenter
     ****************************************/
    
    class Domain1DCenter : public Domain1D {
        
    private:
        
        // Center
        CTptr _center;
        
        // Multiplier
        double _multiplier;
        
        // Internal copy func/clean up
        void _clean_up();
        void _move(Domain1DCenter& other);
        void _copy(const Domain1DCenter& other);
        
    public:
        
        /********************
         Constructor
         ********************/
        
        Domain1DCenter(CTptr center, double multiplier, double delta, double zero);
        Domain1DCenter(const Domain1DCenter& other);
        Domain1DCenter& operator=(const Domain1DCenter& other);
        Domain1DCenter(Domain1DCenter&& other);
        Domain1DCenter& operator=(Domain1DCenter&& other);
        ~Domain1DCenter();
        
        /********************
         Getters
         ********************/
        
        CTptr get_center() const;
        
        double get_val_at_timepoint(int timepoint) const;
    };










    /****************************************
     Domain
     ****************************************/
    
    class Domain {
        
    private:

        // No coeffs
        int _no_dims;
        
        // All my domains
        std::vector<Domain1D*> _domain;
        std::vector<Domain1DParam*> _domain_param; // extra storage but wever
        std::vector<Domain1DCenter*> _domain_center; // extra storage but wever
        std::vector<Domain1DObs*> _domain_obs; // extra storage but wever
        std::map<ITptr, int> _param_deriv_idxs;

        // Vals at different timepoints
        std::map<int,std::vector<double>> _vals;
        std::vector<double> _val_substep;
        
        // Internal copy func/clean up
        void _clean_up();
        void _move(Domain& other);
        void _copy(const Domain& other);
        void _shared_constructor(std::vector<Domain1D*> domain);

    public:
        
        /********************
         Constructor
         ********************/
        
        Domain(std::vector<Domain1DParam*> domains);
        Domain(std::vector<Domain1DCenter*> domains);
        Domain(std::vector<Domain1DObs*> domains);
        Domain(const Domain& other);
        Domain& operator=(const Domain& other);
        Domain(Domain&& other);
        Domain& operator=(Domain&& other);
        ~Domain();
        
        /********************
         Getters
         ********************/
        
        int get_no_dims() const;
        
        const std::vector<Domain1D*>& get_domain() const;
        const std::vector<Domain1DParam*>& get_domain_param() const;
        const std::vector<Domain1DCenter*>& get_domain_center() const;
        const std::vector<Domain1DObs*>& get_domain_obs() const;

        void calculate_val_at_timepoint(int timepoint);
        void calculate_substep_val();

        const std::vector<double>& get_val_at_timepoint(int timepoint) const;
        const std::vector<double>& get_substep_val() const;

        int* get_idx_of_ixn_param(ITptr ixn) const;
    };



















	/****************************************
	DiffEqRHS
	****************************************/

	class DiffEqRHS : public q3c1::Grid {

    protected:
        
        // No coeffs
        int _no_coeffs;

        // Domain - the domain in dcu::Grid does not store the ixn funcs
        std::shared_ptr<Domain> _domain;

        // Parent ixn param
        ITptr _parent_ixn_param_traj;
        
        // Updates
        std::map<q3c1::Vertex*,std::vector<double>> _updates;

	private:

        // Learning rate
        double _lr;
        
		// Name
		std::string _name;

		// No dims
		int _no_dims;
        
        // Coefficient order
        std::vector<std::vector<q3c1::DimType>> _coeff_order;
        
		// Nesterov
		std::map<q3c1::Vertex*,std::vector<double>> *_nesterov_y_s, *_nesterov_y_sp1;
        
        // Adam
        std::map<q3c1::Vertex*,std::vector<double>> *_adam_m, *_adam_v;

        // Maximum magnitude for update to coeffs
        double *_mag_max_update;
        
		// Internal copy/clean up function
		void _copy(const DiffEqRHS& other);
		void _move(DiffEqRHS& other);
		void _clean_up();

	public:

		/********************
		Constructor
		********************/

		// Note: ownership of domain is NOT transferred
        DiffEqRHS(std::string name, ITptr parent_ixn_param_traj, std::shared_ptr<Domain> domain, double lr);
        DiffEqRHS(const DiffEqRHS& other);
		DiffEqRHS(DiffEqRHS&& other);
		DiffEqRHS& operator=(const DiffEqRHS& other);
		DiffEqRHS& operator=(DiffEqRHS &&other);
		virtual ~DiffEqRHS();
        
        // ***************
        // MARK: - Learning rate
        // ***************
        
        void set_lr(double lr);
        double get_lr() const;

		/********************
		Validate setup
		********************/

		void validate_setup() const;

        // ***************
        // MARK: - Set magnitude of max update
        // ***************
        
        void set_mag_max_update(double mag);
        
		/********************
		Getters
		********************/
        
        std::shared_ptr<Domain> get_domain() const;

		std::string get_name() const;

		ITptr get_parent_ixn_param_traj() const;

        // Get cell at timepoint
        q3c1::Cell* get_cell_at_timepoint(int timepoint, bool form_abscissas=true) const;
        
        // Get val
		double get_val_at_timepoint(int timepoint, bool form_abscissas=true) const;
        double get_substep_val(bool form_abscissas=true) const;

		// Deriv wrt specific coefficient of some basis
		double get_deriv_wrt_u_at_timepoint(int timepoint, q3c1::IdxSet global_vertex_idxs, std::vector<q3c1::DimType> dim_types, bool form_abscissas=true) const;

		// Spatial deriv
		double get_deriv_wrt_nu_at_timepoint(int timepoint, int deriv_dim, bool form_abscissas=true) const;
        double get_deriv_wrt_nu_at_timepoint(int timepoint, ITptr deriv_ixn_param, bool form_abscissas=true) const;
        
        // ***************
        // MARK: - Fix vertices at some timepoint
        // ***************
        
        void fix_all_verts_around_at_timepoint(int timepoint, bool fixed, bool form_abscissas=true) const;
        
		/********************
		Update
		********************/

		// Calculate the update
		// t_start = inclusive
		// t_end = inclusive
        virtual void update_calculate_and_store(int timepoint_start, int timepoint_end, double dt, bool form_abscissas=true);

		// Committ the update
        void update_committ_stored_sgd();
        void update_committ_stored_nesterov(double nesterov_acc);
        void update_committ_stored_adam(int opt_step, double beta_1, double beta_2, double eps);

        // Verbose
        void print_update_stored() const;
    };

    
    /****************************************
     DiffEqRHSCenteredHomWeight
     ****************************************/
    
    class Adjoint;
    
    class DiffEqRHSCenteredHomWeight : public DiffEqRHS {
        
    private:

        // Conn mult
        int _conn_mult;
        
        // Biases adjoints
        std::shared_ptr<Adjoint> _bias_lower_adjoint, _bias_upper_adjoint;
        
        // Centers
        CTptr _center_lower, _center_upper;
        
        // Internal copy/clean up function
        void _copy(const DiffEqRHSCenteredHomWeight& other);
        void _move(DiffEqRHSCenteredHomWeight& other);
        void _clean_up();
        
    public:
        
        /********************
         Constructor
         ********************/
        
        // Note: ownership of domain is NOT transferred
        DiffEqRHSCenteredHomWeight(std::string name, ITptr parent_ixn_param_traj, std::shared_ptr<Domain> domain, double lr, int conn_mult, std::shared_ptr<Adjoint> bias_lower_adjoint, std::shared_ptr<Adjoint> bias_upper_adjoint, CTptr center_lower, CTptr center_upper);
        DiffEqRHSCenteredHomWeight(const DiffEqRHSCenteredHomWeight& other);
        DiffEqRHSCenteredHomWeight(DiffEqRHSCenteredHomWeight&& other);
        DiffEqRHSCenteredHomWeight& operator=(const DiffEqRHSCenteredHomWeight& other);
        DiffEqRHSCenteredHomWeight& operator=(DiffEqRHSCenteredHomWeight &&other);
        ~DiffEqRHSCenteredHomWeight();

        /********************
         Update
         ********************/
        
        // Calculate the update
        // t_start = inclusive
        // t_end = inclusive
        void update_calculate_and_store(int timepoint_start, int timepoint_end, double dt, bool form_abscissas=true);
    };
};
