#include <vector>
#include <string>

#ifndef FWDS_SPECIES_H
#define FWDS_SPECIES_H
#include "fwds/fwds_species.hpp"
#endif

#ifndef FWDS_IXN_PARAM_H
#define FWDS_IXN_PARAM_H
#include "fwds/fwds_ixn_param.hpp"
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	// Forwards
	class DiffEqRHS;
	class Moment;
	class Adjoint;

	/****************************************
	Interaction parameter
	****************************************/

	// Enumeration of type of dimension
	enum class IxnParamType: unsigned int { H, J, K, W, B };

	class IxnParam {

	private:

		class Impl;
		std::unique_ptr<Impl> _impl;

	public:

		/********************
		Constructor
		********************/

		IxnParam(std::string name, IxnParamType type, double init_cond); 
		IxnParam(const IxnParam& other);
		IxnParam(IxnParam&& other);
		IxnParam& operator=(const IxnParam& other);
		IxnParam& operator=(IxnParam&& other);
		~IxnParam();
		
		/********************
		Timesteps
		********************/

		int get_no_timesteps() const;
		void set_no_timesteps(int no_timesteps);

		/********************
		Init cond
		********************/

		double get_init_cond() const;
		void set_init_cond(double init_cond);

		/********************
		Fixed value to IC
		********************/

		void set_fix_value_to_init_cond(bool fixed);
		bool get_is_val_fixed() const;

		/********************
		Name, type
		********************/

		std::string get_name() const;

		IxnParamType get_type() const;

		/********************
		Value
		********************/

		double get_val_at_timepoint(int timepoint) const;

		/********************
		Diff eq
		********************/

		std::shared_ptr<DiffEqRHS> get_diff_eq_rhs() const;
		void set_diff_eq_rhs(std::shared_ptr<DiffEqRHS> diff_eq);

		void solve_diff_eq_at_timepoint_to_plus_one(int timepoint, double dt);

		/********************
		Moment
		********************/

		std::shared_ptr<Moment> get_moment() const;

		/********************
		Adjoint
		********************/

		void set_adjoint(std::shared_ptr<Adjoint> adjoint);
		std::shared_ptr<Adjoint> get_adjoint() const;

		/********************
		Write to file
		********************/

		void write_to_file(std::string fname) const;
	};

};

