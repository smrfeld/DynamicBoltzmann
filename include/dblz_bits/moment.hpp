#include <vector>
#include <string>
#include <map>

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

		// No chains
        std::map<MCType, int> _no_markov_chains;

        // Reaped values from the sampler
		// Size = _batch_size
        std::map<MCType, std::vector<double>> _vals_reaped;

		// Averaged values
        std::map<MCType, double> _val_averaged;
        
        // Difference
        double _val_diff;
	
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
		Batch size/no markov chains
		********************/

        int get_no_markov_chains(MCType type) const;
        void set_no_markov_chains(MCType type, int no_markov_chains);
        
		/********************
		Reset
		********************/

		void reset_moment_samples_to_zero();

		/********************
		Fixed awake
		********************/

		void set_is_awake_moment_fixed(bool flag, double val);
		bool get_is_awake_moment_fixed() const;

		/********************
		Get/set moment
		********************/
        
        // Get moment
        double get_moment(MCType type) const;
        
		// Sample
		void set_moment_sample(MCType type, int i_sample, double val);
        // void increment_moment_sample(MCType type, int i_sample, double val);

		// Average reaps
		void average_moment_samples();
        
        // Get moment difference
        double get_moment_diff_awake_minus_asleep() const;
        
        // Augment moment difference by some value
        void increment_moment_diff_awake_minus_asleep(double val);
        
		/********************
		Write
		********************/

		void write_to_file(std::string fname, bool append=false) const;
	};
};
