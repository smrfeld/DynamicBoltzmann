#include "../include/dblz_bits/moment_diff.hpp"

// Other headers
#include "../include/dblz_bits/general.hpp"
#include "../include/dblz_bits/ixn_param.hpp"
#include "../include/dblz_bits/species.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

/************************************
* Namespace for bmla
************************************/

namespace dblz {

	/****************************************
	MomentDiff
	****************************************/
	
	/********************
	Constructor
	********************/

	MomentDiff::MomentDiff(std::string name, IxnParamType type) {
		_name = name;
		_type = type;

		// Fixed awake moment
		_is_awake_moment_fixed = false;
        
        // Offset
        _val_diff_offset = 0.0;
        
        // Vals
        _val_averaged[MCType::AWAKE] = 0.0;
        _val_averaged[MCType::ASLEEP] = 0.0;
    };
	MomentDiff::MomentDiff(const MomentDiff& other) {
		_copy(other);
	};
	MomentDiff::MomentDiff(MomentDiff&& other) {
		_move(other);
	};
	MomentDiff& MomentDiff::operator=(const MomentDiff& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	MomentDiff& MomentDiff::operator=(MomentDiff&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
	};
	MomentDiff::~MomentDiff() {
		_clean_up();
	};

	void MomentDiff::_clean_up() {
    };
	void MomentDiff::_move(MomentDiff &other) {
		_name = other._name;
		_type = other._type;

		// averaged
		_val_averaged = other._val_averaged;
        
        _val_diff_offset = other._val_diff_offset;
        
		_is_awake_moment_fixed = other._is_awake_moment_fixed;

		// Reset the other
		other._name = "";

		other._val_averaged[MCType::AWAKE] = 0.0;
        other._val_averaged[MCType::ASLEEP] = 0.0;
        
        other._val_diff_offset = 0.0;
        
		other._is_awake_moment_fixed = false;
	};
	void MomentDiff::_copy(const MomentDiff& other) {
		_name = other._name;
		_type = other._type;
        
        // _no_markov_chains = other._no_markov_chains;
        
		// reaped
        // _vals_reaped = other._vals_reaped;

		// averaged
		_val_averaged = other._val_averaged;
        
        _val_diff_offset = other._val_diff_offset;

		_is_awake_moment_fixed = other._is_awake_moment_fixed;
	};

	/********************
	Verbose
	********************/

	void MomentDiff::print_moment_comparison() const {
        std::cout << "(" << _val_averaged.at(MCType::AWAKE) << "," << _val_averaged.at(MCType::ASLEEP) << ") " << std::endl;
	};
    
    std::string MomentDiff::get_moment_comparison_str() const {
        std::stringstream s;
        s << "(" << _val_averaged.at(MCType::AWAKE) << "," << _val_averaged.at(MCType::ASLEEP) << ")";
        return s.str();
    };
    
	/********************
	Name
	********************/

	std::string MomentDiff::get_name() const {
		return _name;
	};
	IxnParamType MomentDiff::get_type() const {
		return _type;
	};

	/********************
	Fixed awake
	********************/

	void MomentDiff::set_is_awake_moment_fixed(bool flag, double val) {
		_is_awake_moment_fixed = flag;
        
        _val_averaged[MCType::AWAKE] = val;
	};
	bool MomentDiff::get_is_awake_moment_fixed() const {
		return _is_awake_moment_fixed;
	};
    
	/********************
	Get/set moment
	********************/

    // Get moment
    double MomentDiff::get_moment(MCType type) const {
        return _val_averaged.at(type);
    };
    void MomentDiff::increment_moment(MCType type, double val) {
        _val_averaged[type] += val;
    };
    void MomentDiff::set_moment(MCType type, double val) {
        _val_averaged[type] = val;
    };
    void MomentDiff::reset_moment(MCType type) {
        _val_averaged[type] = 0.0;
    };

    // Augment moment difference by some value
    void MomentDiff::set_moment_offset(double val) {
        _val_diff_offset = val;
    };

    // Get moment difference
    double MomentDiff::get_moment_diff_awake_minus_asleep_plus_offset() const {
        return _val_averaged.at(MCType::AWAKE) - _val_averaged.at(MCType::ASLEEP) + _val_diff_offset;
    };
    
	/********************
	Write
	********************/

	void MomentDiff::write_to_file(std::string fname, bool append) const {
		std::ofstream f;

		// Open
		if (append) {
			f.open(fname, std::ofstream::out | std::ofstream::app);
		} else {
			f.open(fname);
		};

		// Make sure we found it
		if (!f.is_open()) {
			std::cerr << ">>> Error: MomentDiff::write_to_file <<< could not write to file: " << fname << std::endl;
			exit(EXIT_FAILURE);
		};

        f << _val_averaged.at(MCType::AWAKE) << " " << _val_averaged.at(MCType::ASLEEP) << "\n";

		// Close
		f.close();
    };
};


