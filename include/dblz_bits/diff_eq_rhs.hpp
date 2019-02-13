#include <vector>
#include <string>
#include <map>

#include "fwds/fwds_ixn_param_traj.hpp"

#include <q3c1> // guard is built-in

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	Domain1D
	****************************************/

	class Domain1D : public q3c1::Dimension1D {

	private:

		ITptr _ixn_param_traj;

		// Internal copy func/clean up
		void _clean_up();
		void _move(Domain1D& other);
		void _copy(const Domain1D& other);

	public:

		/********************
		Constructor
		********************/

		Domain1D(ITptr ixn_param_traj, double delta, double zero);
		Domain1D(const Domain1D& other);
		Domain1D& operator=(const Domain1D& other);
		Domain1D(Domain1D&& other);
		Domain1D& operator=(Domain1D&& other);
		~Domain1D();

		/********************
		Getters
		********************/

		std::string get_name() const;
		ITptr get_ixn_param_traj() const;
	};































	/****************************************
	DiffEqRHS
	****************************************/

	class DiffEqRHS : public q3c1::Grid {

	private:

        // Learning rate
        double _lr;
        
		// Name
		std::string _name;

		// No dims
		int _no_dims;

        // No coeffs
        int _no_coeffs;
        
        // Coefficient order
        std::vector<std::vector<q3c1::DimType>> _coeff_order;
        
		// Domain - the domain in dcu::Grid does not store the ixn funcs
		std::vector<Domain1D*> _domain;

		// Parent ixn param
		ITptr _parent_ixn_param_traj;

		// Updates
		std::map<q3c1::Vertex*,std::vector<double>> _updates;

		// Nesterov
		std::map<q3c1::Vertex*,std::vector<double>> *_nesterov_y_s, *_nesterov_y_sp1;
        
        // Adam
        std::map<q3c1::Vertex*,std::vector<double>> *_adam_m, *_adam_v;

		// Internal
        mutable std::vector<double> _abscissas;
        void _form_abscissas(int timepoint) const;

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
		DiffEqRHS(std::string name, ITptr parent_ixn_param_traj, std::vector<Domain1D*> domain, double lr);
		DiffEqRHS(const DiffEqRHS& other);
		DiffEqRHS(DiffEqRHS&& other);
		DiffEqRHS& operator=(const DiffEqRHS& other);
		DiffEqRHS& operator=(DiffEqRHS &&other);
		~DiffEqRHS();
        
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

		std::string get_name() const;

		ITptr get_parent_ixn_param_traj() const;

		const std::vector<Domain1D*>& get_domain() const;

		double get_val_at_timepoint(int timepoint) const;

		// Deriv wrt specific coefficient of some basis
		double get_deriv_wrt_u_at_timepoint(int timepoint, q3c1::IdxSet global_vertex_idxs, std::vector<q3c1::DimType> dim_types) const;

		// Spatial deriv
		double get_deriv_wrt_nu_at_timepoint(int timepoint, int deriv_dim) const;

		/********************
		Update
		********************/

		// Calculate the update
		// t_start = inclusive
		// t_end = inclusive
		// void update_calculate_and_store(int timepoint_start, int timepoint_end, double dt, double dopt);
        void update_calculate_and_store(int timepoint_start, int timepoint_end, double dt);

		// Committ the update
		// void update_committ_stored(bool nesterov_mode=true, double nesterov_acc=0.5);
        void update_committ_stored_sgd();
        void update_committ_stored_nesterov(double nesterov_acc);
        void update_committ_stored_adam(int opt_step, double beta_1, double beta_2, double eps);

        // Verbose
        void print_update_stored() const;
    };

};
