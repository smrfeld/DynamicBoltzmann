#ifndef VECTOR_h
#define VECTOR_h
#include <vector>
#endif

#ifndef STRING_h
#define STRING_h
#include <string>
#endif

#ifndef IXN_PARAM_TRAJ_h
#define IXN_PARAM_TRAJ_h
#include "ixn_param_traj.hpp"
#endif

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {	

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
		std::vector<IxnParamTraj*> _ixn_params;

		// Dimension length squares
		std::vector<int> _dim_pwrs;

	public:

		// Constructor
		Array(std::vector<IxnParamTraj*> ixn_params);
		Array(IxnParamTraj* ixn_param);
		Array(const Array& other);
		Array& operator=(const Array& other);
		Array(Array&& other);
		Array& operator=(Array&& other);
		~Array();

		// Get/set an element by index
		double get_by_idxs(int *idxs) const;
		double get_by_idx(int i) const;
		void set_by_idxs(int *idxs, double val);
		void set_by_idx(int i, double val);

		// Get indexes by element
		void get_idxs(int i, int* idxs) const;

		// Write/Read to a file
		void write_grid(std::string fname) const;
		void write_vals(std::string dir, std::string name, int idx) const;
		void read_vals(std::string fname);

		// Check dimensions against another array
		bool check_dims(const Array& other) const;

		// Zero
		void zero();
	};

	/****************************************
	Forward declare
	****************************************/

	class VarTermTraj;

	/****************************************
	BasisFunc
	****************************************/

	class BasisFunc : public Array {

	private:

		// Name
		std::string _name;

		// The ixn param that I am the basis function for
		IxnParamTraj* _parent_ixn_param;

		// The variational terms and ixn_params needed to update
		std::vector<std::pair<IxnParamTraj*,VarTermTraj*>> _update_ptrs;

		// Get bounding n-dim cube of 4 pts
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
		void _copy(const BasisFunc& other);
		void _clean_up();

		// Interpolation
		double _cubic_interp(double p[4], double f);
		double _deriv(double p[4], double f);
		double _n_cubic_interp(int dim, double* p, double fracs[], bool derivs[]);

	public:

		// Constructor
		BasisFunc(std::string name, IxnParamTraj* parent_ixn_param, std::vector<IxnParamTraj*> ixn_params);
		BasisFunc(const BasisFunc& other);
		BasisFunc& operator=(const BasisFunc& other);
		~BasisFunc();

		// Move to the nesterov intermediate point
		void nesterov_move_to_intermediate_pt(int opt_step);

		// Set prev nesterov
		void nesterov_set_prev_equal_curr();

		// Add pointers needed to update
		void add_update_ptrs(IxnParamTraj* ixn_param, VarTermTraj* var_term);

		// Validate setup
		void validate_setup() const;

		// Get values, if they are in the lattice
		double get_at_time(int it);
		double get_deriv_at_time(int it, int i_dim);

		// Name
		std::string name() const;

		// Calculate the new basis function
		// t_start (inclusive) to end (non-inclusive)
		void update(int t_start, int t_end, double dt, double dopt, bool exp_decay, double exp_decay_t0, double exp_decay_lambda, bool l2_reg_params_mode, std::map<IxnParamTraj*,double> &l2_lambda_params, std::map<IxnParamTraj*,double> &l2_reg_centers); 
		void update_gather(int t_start, int t_end, double dt, double dopt, bool exp_decay, double exp_decay_t0, double exp_decay_lambda);
		void update_committ_gathered(int t_start, int t_end, double dt, double dopt, bool l2_reg_params_mode, std::map<IxnParamTraj*,double> &l2_lambda_params, std::map<IxnParamTraj*,double> &l2_reg_centers);

		// Test fill in various dimensions
		void test_fill_2d();
		void test_fill_3d();

		// From parent
		
		// Get/set an element by index
		double get_by_idxs(int *idxs) const;
		double get_by_idx(int i) const;
		void set_by_idxs(int *idxs, double val);

		// Write/Read grid/vals
		void write_grid(std::string fname) const;
		void write_vals(std::string dir, int idx_opt_step) const;
		void read_vals(std::string fname);

		// Get the delta source
		double get_delta_source(int it, int i);
	};

};