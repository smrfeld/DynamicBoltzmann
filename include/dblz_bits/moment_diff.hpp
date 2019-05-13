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
	MomentDiff
	****************************************/
	
	class MomentDiff {

	private:

		// Name
		std::string _name;

		// Type = H, J, K, B, W
		IxnParamType _type;
        
		// Averaged values
        std::map<MCType, double> _val_averaged;
        
        // Offset to the difference
        double _val_diff_offset;
	
		// If the awake moment is fixed
		bool _is_awake_moment_fixed;

		// Constructor helpers
		void _clean_up();
		void _move(MomentDiff &other);
		void _copy(const MomentDiff& other);

	public:

		/********************
		Constructor
		********************/

		MomentDiff(std::string name, IxnParamType type);
		MomentDiff(const MomentDiff& other);
		MomentDiff(MomentDiff&& other);
		MomentDiff& operator=(const MomentDiff& other);
		MomentDiff& operator=(MomentDiff&& other);
		~MomentDiff();

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

		/********************
		Get/set moment
		********************/
        
        // Get moment
        double get_moment(MCType type) const;
        void increment_moment(MCType type, double val);
        void set_moment(MCType type, double val);
        void reset_moment(MCType type);
        
        // Augment moment difference by some value
        void set_moment_offset(double val);

        // Get moment difference
        double get_moment_diff_awake_minus_asleep() const;
        double get_moment_diff_awake_minus_asleep_plus_offset() const;

		/********************
		Write
		********************/

		void write_to_file(std::string fname, bool append=false) const;
    };
};
