#ifndef STRING_H
#define STRING_H
#include <string>
#endif

#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

/************************************
* Namespace for dboltz
************************************/

namespace dboltz {

	/****************************************
	Variational Term Trajectory
	****************************************/

	// Forward
	class IxnParamTraj;
	class BasisFunc;
	class Array;

	class VarTermTraj {
	private:

		// Name
		std::string _name;

		// Length of trajectory
		int _n_t;

		// Numerator and denominator
		IxnParamTraj *_num;
		BasisFunc *_denom;

		// Whether or not this term deserves a delta source
		bool _delta_source;

		// Basis func corresponding to num
		BasisFunc *_num_bf;

		// Derivs of the numerators basis function necessary for updating
		int _n_ixn_params_in_num_bf;
		double *_num_bf_derivs;

		// Values over time
		std::vector<Array> _vals;

		// Val length of each array
		int _val_len;

		// Pointers to other ixn params and variational terms needed to update this term
		// The interaction parameters that are arguments to _num_bf, so we can take derivatives
		std::vector<VarTermTraj*> _update_var_terms;

		// Clean up/copy
		void _copy(const VarTermTraj& other);
		void _clean_up();

	public:

		// Constructor
		VarTermTraj(std::string name, IxnParamTraj *num, BasisFunc *denom, std::vector<IxnParamTraj*> denom_ixn_params, BasisFunc *num_bf, int n_ixn_params_in_num_bf, int n_t);
		VarTermTraj(const VarTermTraj& other);
		VarTermTraj & operator=(const VarTermTraj& other);
		~VarTermTraj();

		// Set the pointers needed to update this term
		void add_update_ptr(VarTermTraj* var_term);

		// Validate setup
		void validate_setup() const;

		// Get numerator/denominator
		IxnParamTraj* get_numerator_ixn_param_traj() const;
		BasisFunc* get_denominator_basis_func() const;

		// Calculate next timestep
		void calculate_at_time(int it_next, double dt);

		// Getters/setters
		double get_at_time_by_idx(int it, int i);
		std::string name();

		// Write
		void write_vals(std::string dir, int idx_opt_step) const;
	};

};