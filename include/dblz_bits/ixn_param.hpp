#include <vector>
#include <string>

#include "fwds/fwds_species.hpp"
#include "fwds/fwds_ixn_param.hpp"

/************************************
* Namespace for bmla
************************************/

namespace dblz {

	// Forwards
	class Moment;

	/****************************************
	Interaction parameter
	****************************************/

	// Enumeration of type of dimension
	enum class IxnParamType: unsigned int { H, J, K, W, B, X, GAMMA, BETA };

	class IxnParam {

	private:
        
        // Learning rate
        double _lr;
        
        // Values
        double _val;
        double _update;
        
        // nesterov
        double *_nesterov_y_s, *_nesterov_y_sp1;
        
        // adam
        double *_adam_m, *_adam_v;
        
        // Fixed value
        bool _is_val_fixed;
        
        // Moment
        std::shared_ptr<Moment> _moment;
        
        // Copy, clean up
        void _clean_up();
        void _copy(const IxnParam& other);
        void _move(IxnParam &other);

	public:

		/********************
		Constructor
		********************/

		IxnParam(std::string name, IxnParamType type, double init_guess, double lr);
		IxnParam(const IxnParam& other);
		IxnParam(IxnParam&& other);
		IxnParam& operator=(const IxnParam& other);
		IxnParam& operator=(IxnParam&& other);
		~IxnParam();

        // ***************
        // MARK: - Learning rate
        // ***************
        
        double get_lr() const;
        void set_lr(double lr);
        
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

        void update_calculate_and_store_l2(double l2_lambda, double l2_center);
		void update_calculate_and_store();
        
		void update_committ_stored_sgd();
		void update_committ_stored_nesterov(double nesterov_acc);
		void update_committ_stored_adam(int opt_step, double beta_1, double beta_2, double eps);

		/********************
		Write to file
		********************/

		void write_to_file(std::string fname, bool append=false) const;
	};

};

