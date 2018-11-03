#include <vector>
#include <string>
#include <map>

#ifndef FWDS_IXN_PARAM_H
#define FWDS_IXN_PARAM_H
#include "fwds/fwds_ixn_param.hpp"
#endif

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
	DiffEqRHS
	****************************************/

	class DiffEqRHS : public q3c1::Grid {

	private:

		// Name
		std::string _name;

		// No dims
		int _no_dims;

		// Domain - the domain in dcu::Grid does not store the ixn funcs
		std::vector<Domain1D*> _domain;

		// Parent ixn param
		Iptr _parent_ixn_param;

		// Updates
		std::map<q3c1::Vertex*,std::vector<double>> _updates;

		// Internal
		std::vector<double> _form_abscissas(int timepoint) const;

		// Internal copy/clean up function
		void _copy(const DiffEqRHS& other);
		void _move(DiffEqRHS& other);
		void _clean_up();

	public:

		/********************
		Constructor
		********************/

		// Note: ownership of domain is NOT transferred
		DiffEqRHS(std::string name, Iptr parent_ixn_param, std::vector<Domain1D*> domain);
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

		const std::vector<Domain1D*>& get_domain() const;

		bool check_val_is_in_domain_at_timepoint(int timepoint) const;

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
		// t_end = non-inclusive
		void update_calculate_and_store(int timepoint_start, int timepoint_end, double dt, double dopt);

		// Verbose
		void print_update_stored() const;

		// Committ the update
		void update_committ_stored(int opt_step, bool nesterov_mode=true);
	};

};