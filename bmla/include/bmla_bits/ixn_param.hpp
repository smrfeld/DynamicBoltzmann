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
* Namespace for bmla
************************************/

namespace bmla {

	// Forwards
	class Moment;

	/****************************************
	Interaction parameter
	****************************************/

	// Enumeration of type of dimension
	enum class IxnParamType: unsigned int { H, J, K, W, B, X };

	class IxnParam {

	private:

		class Impl;
		std::unique_ptr<Impl> _impl;

	public:

		/********************
		Constructor
		********************/

		IxnParam(std::string name, IxnParamType type, double init_guess); 
		IxnParam(const IxnParam& other);
		IxnParam(IxnParam&& other);
		IxnParam& operator=(const IxnParam& other);
		IxnParam& operator=(IxnParam&& other);
		~IxnParam();

		/********************
		Name, type
		********************/

		std::string get_name() const;

		IxnParamType get_type() const;

		/********************
		Value
		********************/

		double get_val() const;
		void set_val(double val);

		/********************
		Fixed value
		********************/

		void set_fix_value(bool fixed);
		bool get_is_val_fixed() const;

		/********************
		Moment
		********************/

		std::shared_ptr<Moment> get_moment() const;

		/********************
		Update
		********************/

		void update_calculate_and_store(double dopt, bool l2_mode=false, double l2_lambda=0.0, double l2_center=0.0);
		void update_committ_stored(bool nesterov_mode=true);

		/********************
		Write to file
		********************/

		void write_to_file(std::string fname, bool append=false) const;
	};

};

