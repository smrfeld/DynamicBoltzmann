#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

#ifndef STRING_H
#define STRING_H
#include <string>
#endif

#ifndef MAP_H
#define MAP_H
#include <map>
#endif

#ifndef FWDS_IXN_PARAM_H
#define FWDS_IXN_PARAM_H
#include "fwds/fwds_ixn_param.hpp"
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	Domain1D
	****************************************/

	class Domain1D {

	private:

		Iptr _ixn_param;
		double _min;
		double _max;
		double _delta;
		int _no_pts;

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
		int get_no_pts() const;
		double get_min() const;
		double get_max() const;
		double get_delta() const;

		// Get pt in domain
		double get_pt_by_idx(int i) const;

		// Check if point is in domain
		bool check_if_pt_is_inside_domain(double x) const;

		// Get indexes surrounding a point
		// ie point is between i and i+1 where i is returned
		int get_idxs_surrounding_pt(double x) const; 

		// Get fraction of a point between successive points
		double get_frac_between(double x) const;
		// Second optional specification: the return of the surrounding idxs
		double get_frac_between(double x, int i) const;

		// Print domain range
		void print_bounds() const;
	};
























	/****************************************
	Array
	****************************************/

	class Array {
	private:
		
		// Internal copy func/clean up
		void _copy(const Array& other);
		void _clean_up();

	protected:

		// Number of dimensions
		int _n_params;

		// Vals
		double *_vals;
		int _val_len;

		// Dimensions
		std::vector<Domain1D> _domain;

		// Dimension length squares
		std::vector<int> _dim_pwrs;

	public:

		/********************
		Constructor
		********************/

		Array(std::vector<Domain1D> domain);
		Array(const Array& other);
		Array& operator=(const Array& other);
		Array(Array&& other);
		Array& operator=(Array&& other);
		~Array();

		/********************
		Get/set an element by index
		********************/

		double get_by_idxs(int *idxs) const;
		double get_by_idx(int i) const;
		void set_by_idxs(int *idxs, double val);
		void set_by_idx(int i, double val);

		/********************
		Get indexes by element
		********************/

		void get_idxs(int i, int* idxs) const;

		/********************
		Write/Read to a file
		********************/

		void write_grid(std::string fname) const;
		void write_vals(std::string fname) const;
		void read_vals(std::string fname);

		/********************
		Check dims against another array
		********************/

		bool check_dims(const Array& other) const;

		/********************
		Reset to zero
		********************/

		void reset_to_zero();
	};




































	/****************************************
	DiffEqRHS
	****************************************/

	class DiffEqRHS : public Array {

	private:

		// Name
		std::string _name;

		// Get bounding n-dim cube (in 2D, this means 4 pts, etc)
		void _get_bounding(int it, bool safe=false);
		void _fill_p(int dim);
		// Parameters - only alloc/dealloc only once
		int* _idxs_bounding;
		int* _idxs_p_cube;
		int* _idxs_ext_1;
		int* _idxs_ext_2;
		double *_fracs;
		double *_p_cube;

		// Update, if needed
		double *_update_gathered;

		// Stored update for nesterov
		Array *_nesterov_prev_pt;

		// Store derivatives in each direction

		// Specify in which dimension to take a derivative
		bool *_derivs;

		// Internal copy/clean up function
		void _copy(const DiffEqRHS& other);
		void _clean_up();

		// Interpolation
		double _cubic_interp(double p[4], double f);
		double _deriv(double p[4], double f);
		double _n_cubic_interp(int dim, double* p, double fracs[], bool derivs[]);

	public:

		/********************
		Constructor
		********************/

		DiffEqRHS(std::string name, std::vector<Domain1D> domain);
		DiffEqRHS(const DiffEqRHS& other);
		DiffEqRHS& operator=(const DiffEqRHS& other);
		~DiffEqRHS();

		/********************
		Validate setup
		********************/

		void validate_setup() const;

		/********************
		Getters
		********************/

		std::string get_name() const;

		double get_val_at_timepoint(int it);
		double get_deriv_at_timepoint(int it, int i_dim);

		/********************
		Nesterov
		********************/

		// Move to the nesterov intermediate point
		void nesterov_move_to_intermediate_pt(int opt_step);

		// Set prev nesterov
		void nesterov_set_prev_equal_curr();

		/********************
		Update
		********************/

		// Calculate the new basis function
		// t_start (inclusive) to end (non-inclusive)
		void update(int t_start, int t_end, double dt, double dopt, bool exp_decay, double exp_decay_t0, double exp_decay_lambda, bool l2_reg_params_mode, std::map<IxnParam*,double> &l2_lambda_params, std::map<IxnParam*,double> &l2_reg_centers); 
		void update_gather(int t_start, int t_end, double dt, double dopt, bool exp_decay, double exp_decay_t0, double exp_decay_lambda, bool l2_reg_params_mode, std::map<IxnParam*,double> &l2_lambda_params, std::map<IxnParam*,double> &l2_reg_centers);
		void update_committ_gathered();

	};

};