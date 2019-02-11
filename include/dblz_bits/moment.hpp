#include <vector>
#include <string>
#include <map>
#include <armadillo>

#include "fwds/fwds_species.hpp"

/************************************
* Namespace for bmla
************************************/

namespace dblz {

	// Forward
    enum class IxnParamType: unsigned int;
    enum class MCType: unsigned int { AWAKE, ASLEEP };

	/****************************************
	Moment
	****************************************/
	
	class Moment {

	private:

		// Name
		std::string _name;

		// Type = H, J, K, B, W
		IxnParamType _type;
        
		// Averaged values
        std::map<MCType, double> _val_averaged;
        
        // Average mat
        // Only for W or X
        std::map<MCType,arma::mat*> _weight_matrix;
        arma::mat* _weight_matrix_awake_minus_asleep;
        // Only for biases
        std::map<MCType,arma::vec*> _bias_vec;
        arma::vec* _bias_vec_awake_minus_asleep;

        // Offset to the difference
        double _val_diff_offset;
	
		// If the awake moment is fixed
		bool _is_awake_moment_fixed;

		// Constructor helpers
		void _clean_up();
		void _move(Moment &other);
		void _copy(const Moment& other);

	public:

		/********************
		Constructor
		********************/

		Moment(std::string name, IxnParamType type);
		Moment(const Moment& other);
		Moment(Moment&& other);
		Moment& operator=(const Moment& other);
		Moment& operator=(Moment&& other);
		~Moment();

		/********************
		Verbose
		********************/

		void print_moment_comparison() const;
        std::string get_moment_comparison_str() const;

		/********************
		Name, type
		********************/

		std::string get_name() const;
		IxnParamType get_type() const;

		/********************
		Fixed awake
		********************/

		void set_is_awake_moment_fixed(bool flag, double val);
		bool get_is_awake_moment_fixed() const;

        // ***************
        // MARK: - Set W matrix / bias vec
        // ***************
        
        // Weight matrix is always lower layer to higher layer
        void set_weight_matrix_dims(int dim_upper_layer, int dim_lower_layer);
        void reset_weight_matrix(MCType type);
        void increment_weight_matrix(MCType type, arma::mat increment);
        void set_moment_to_weight_matrix_sum(MCType type);
        const arma::mat& get_weight_matrix(MCType type) const;

        void calculate_weight_matrix_awake_minus_asleep();
        const arma::mat& get_weight_matrix_awake_minus_asleep() const;
        
        void set_bias_vec_dims(int dims);
        void reset_bias_vec(MCType type);
        void increment_bias_vec(MCType type, arma::vec increment);
        void set_moment_to_bias_vec_sum(MCType type);
        const arma::vec& get_bias_vec(MCType type) const;

        void calculate_bias_vec_awake_minus_asleep();
        const arma::vec& get_bias_vec_awake_minus_asleep() const;

		/********************
		Get/set moment
		********************/
        
        // Get moment
        double get_moment(MCType type) const;
        void increment_moment(MCType type, double val);
        void set_moment(MCType type, double val);
        void reset_moment(MCType type);
        
        // Get moment difference
        double get_moment_diff_awake_minus_asleep() const;
        
        // Augment moment difference by some value
        void set_moment_diff_awake_minus_asleep_offset(double val);
        
		/********************
		Write
		********************/

		void write_to_file(std::string fname, bool append=false) const;
	};
};
