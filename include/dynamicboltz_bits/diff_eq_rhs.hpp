#include <vector>
#include <string>
#include <map>

#ifndef FWDS_IXN_PARAM_H
#define FWDS_IXN_PARAM_H
#include "fwds/fwds_ixn_param.hpp"
#endif

#include <dcubic> // guard is built-in

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	Domain1D
	****************************************/

	class Domain1D : public dcu::Dimension1D {

	private:

		Iptr _ixn_param;

		// Internal copy func/clean up
		void _clean_up();
		void _reset();
		void _copy(const Domain1D& other);

	public:

		/********************
		Constructor
		********************/

		Domain1D(Iptr ixn_param, double min, double max, int no_pts);
		Domain1D(const Domain1D& other);
		Domain1D& operator=(const Domain1D& other);
		Domain1D(Domain1D&& other);
		Domain1D& operator=(Domain1D&& other);
		~Domain1D();

		/********************
		Getters
		********************/

		std::string get_name() const;
		Iptr get_ixn_param() const;
	};


























	/****************************************
	Domain
	****************************************/

	class Domain {

	private:

		std::vector<dcu::Dimension1D> _dimensions;
		std::vector<Domain1D> _domain;

		// Internal copy func/clean up
		void _clean_up();
		void _reset();
		void _copy(const Domain& other);

	public:

		/********************
		Constructor
		********************/

		Domain();
		Domain(std::vector<Domain1D> domain);
		Domain(const Domain& other);
		Domain& operator=(const Domain& other);
		Domain(Domain&& other);
		Domain& operator=(Domain&& other);
		~Domain();

		/********************
		Setters
		********************/

		void add_dimension(Domain1D domain);

		/********************
		Getters
		********************/

		int size() const;
		std::vector<dcu::Dimension1D> get_dimensions() const;
		std::vector<Domain1D> get_domain() const;
	};































	/****************************************
	DiffEqRHS
	****************************************/

	class DiffEqRHS : public dcu::Grid {

	private:

		// Name
		std::string _name;

		// Domain - the domain in dcu::Grid does not store the ixn funcs
		std::vector<Domain1D> _domain;
		std::vector<dcu::Dimension1D> _dimensions;

		// Helper structures for evaluating
		std::vector<double> _abscissas;
		dcu::IdxSet4 _idxs_k;

		// Parent ixn param
		Iptr _parent_ixn_param;

		// Updates
		std::map<dcu::GridPtKey, double> _updates;
		dcu::Grid *_nesterov_prev_pt;

		// Internal
		void _form_abscissas(int timepoint);

		// Internal copy/clean up function
		void _copy(const DiffEqRHS& other);
		void _move(DiffEqRHS& other);
		void _clean_up();

	public:

		/********************
		Constructor
		********************/

		DiffEqRHS(std::string name, Iptr parent_ixn_param, Domain domain);
		DiffEqRHS(const DiffEqRHS& other);
		DiffEqRHS(DiffEqRHS&& other);
		DiffEqRHS& operator=(const DiffEqRHS& other);
		DiffEqRHS& operator=(DiffEqRHS &&other);
		~DiffEqRHS();

		/********************
		Validate setup
		********************/

		void validate_setup() const;

		/********************
		Getters
		********************/

		std::string get_name() const;

		Iptr get_parent_ixn_param() const;

		const std::vector<Domain1D>& get_domain() const;

		double get_val_at_timepoint(int timepoint);
		double get_deriv_wrt_u_at_timepoint(int timepoint, dcu::IdxSet4 idxs_k);
		double get_deriv_wrt_nu_at_timepoint(int timepoint, int k);

		/********************
		Nesterov
		********************/

		// Move to the nesterov intermediate point
		void nesterov_move_to_intermediate_pt(int i_opt_step);

		// Set prev nesterov
		void nesterov_set_prev_equal_curr();

		/********************
		Update
		********************/

		// Calculate the update
		// t_start = inclusive
		// t_end = non-inclusive
		void update_calculate_and_store(int timepoint_start, int timepoint_end, double dt, double dopt, bool clear_existing_updates=true);

		// Verbose
		void print_update_stored() const;

		// Committ the update
		void update_committ_stored();

		// Calculate the new basis function
		// t_start (inclusive) to end (non-inclusive)
		/*
		void update(int t_start, int t_end, double dt, double dopt, bool exp_decay, double exp_decay_t0, double exp_decay_lambda, bool l2_reg_params_mode, std::map<IxnParam*,double> &l2_lambda_params, std::map<IxnParam*,double> &l2_reg_centers); 
		void update_gather(int t_start, int t_end, double dt, double dopt, bool exp_decay, double exp_decay_t0, double exp_decay_lambda, bool l2_reg_params_mode, std::map<IxnParam*,double> &l2_lambda_params, std::map<IxnParam*,double> &l2_reg_centers);
		void update_committ_gathered();
		*/
	};

};