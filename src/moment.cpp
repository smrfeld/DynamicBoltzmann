#include "../include/dblz_bits/moment.hpp"

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
	Moment
	****************************************/
	
	/********************
	Constructor
	********************/

	Moment::Moment(std::string name, IxnParamType type) {
		_name = name;
		_type = type;

        set_no_markov_chains(MCType::AWAKE,1);
        set_no_markov_chains(MCType::ASLEEP,1);

		// Fixed awake moment
		_is_awake_moment_fixed = false;
        
        // Offset
        _val_diff_offset = 0.0;
	};
	Moment::Moment(const Moment& other) {
		_copy(other);
	};
	Moment::Moment(Moment&& other) {
		_move(other);
	};
	Moment& Moment::operator=(const Moment& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	Moment& Moment::operator=(Moment&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
	};
	Moment::~Moment() {
		_clean_up();
	};

	void Moment::_clean_up() {};
	void Moment::_move(Moment &other) {
		_name = other._name;
		_type = other._type;

        _no_markov_chains = other._no_markov_chains;
        
		// reaped
		_vals_reaped = other._vals_reaped;

		// averaged
		_val_averaged = other._val_averaged;
		
        _val_diff_offset = other._val_diff_offset;
        
		_is_awake_moment_fixed = other._is_awake_moment_fixed;

		// Reset the other
		other._name = "";

        other._no_markov_chains[MCType::AWAKE] = 0;
        other._no_markov_chains[MCType::ASLEEP] = 0;

		other._vals_reaped[MCType::AWAKE].clear();
        other._vals_reaped[MCType::ASLEEP].clear();

		other._val_averaged[MCType::AWAKE] = 0.0;
        other._val_averaged[MCType::ASLEEP] = 0.0;

        other._val_diff_offset = 0.0;
        
		other._is_awake_moment_fixed = false;
	};
	void Moment::_copy(const Moment& other) {
		_name = other._name;
		_type = other._type;
        
        _no_markov_chains = other._no_markov_chains;
        
		// reaped
        _vals_reaped = other._vals_reaped;

		// averaged
		_val_averaged = other._val_averaged;
        
        _val_diff_offset = other._val_diff_offset;

		_is_awake_moment_fixed = other._is_awake_moment_fixed;
	};

	/********************
	Verbose
	********************/

	void Moment::print_moment_comparison() const {
        std::cout << "(" << _val_averaged.at(MCType::AWAKE) << "," << _val_averaged.at(MCType::ASLEEP) << ") " << std::endl;
	};
    
    std::string Moment::get_moment_comparison_str() const {
        std::stringstream s;
        s << "(" << _val_averaged.at(MCType::AWAKE) << "," << _val_averaged.at(MCType::ASLEEP) << ")";
        return s.str();
    };
    
	/********************
	Name
	********************/

	std::string Moment::get_name() const {
		return _name;
	};
	IxnParamType Moment::get_type() const {
		return _type;
	};

	/********************
	Batch size
	********************/

    int Moment::get_no_markov_chains(MCType type) const {
        return _no_markov_chains.at(type);
    };
    void Moment::set_no_markov_chains(MCType type, int no_markov_chains) {
        // Params
        _no_markov_chains[type] = no_markov_chains;
        
        // Reaped vals
        _vals_reaped[type] = std::vector<double>(_no_markov_chains[type]);
        
        // Reset averaged vals
        _val_averaged[type] = 0.0;
    };

    
	/********************
	Reset
	********************/

	void Moment::reset_moment_samples_to_zero() {
        for (auto i_chain=0; i_chain<_no_markov_chains[MCType::AWAKE]; i_chain++) {
            set_moment_sample(MCType::AWAKE, i_chain, 0.0);
        };
        for (auto i_chain=0; i_chain<_no_markov_chains[MCType::ASLEEP]; i_chain++) {
            set_moment_sample(MCType::ASLEEP, i_chain, 0.0);
        };
	};

	/********************
	Fixed awake
	********************/

	void Moment::set_is_awake_moment_fixed(bool flag, double val) {
		_is_awake_moment_fixed = flag;
        
        _val_averaged[MCType::AWAKE] = val;
	};
	bool Moment::get_is_awake_moment_fixed() const {
		return _is_awake_moment_fixed;
	};

	/********************
	Get/set moment
	********************/

    // Get moment
    double Moment::get_moment(MCType type) const {
        return _val_averaged.at(type);
    };
    
	// Batch
	void Moment::set_moment_sample(MCType type, int i_sample, double val) {
        _vals_reaped[type][i_sample] = val;
	};
    /*
    void Moment::increment_moment_sample(MCType type, int i_sample, double val) {
        _vals_reaped[type][i_sample] += val;
    };
     */

    // Average reaps
	void Moment::average_moment_samples() {
        _val_averaged[MCType::AWAKE] = 0.0;
        _val_averaged[MCType::ASLEEP] = 0.0;
        for (auto i=0; i<_no_markov_chains[MCType::AWAKE]; i++) {
            _val_averaged[MCType::AWAKE] += _vals_reaped[MCType::AWAKE][i];
        };
        for (auto i=0; i<_no_markov_chains[MCType::ASLEEP]; i++) {
            _val_averaged[MCType::ASLEEP] += _vals_reaped[MCType::ASLEEP][i];
        };
        _val_averaged[MCType::AWAKE] /= _no_markov_chains[MCType::AWAKE];
        _val_averaged[MCType::ASLEEP] /= _no_markov_chains[MCType::ASLEEP];
    };

    // Get moment difference
    double Moment::get_moment_diff_awake_minus_asleep() const {
        return _val_averaged.at(MCType::AWAKE) - _val_averaged.at(MCType::ASLEEP) + _val_diff_offset;
    };
    
    // Augment moment difference by some value
    void Moment::set_moment_diff_awake_minus_asleep_offset(double val) {
        _val_diff_offset = val;
    };
    
	/********************
	Write
	********************/

	void Moment::write_to_file(std::string fname, bool append) const {
		std::ofstream f;

		// Open
		if (append) {
			f.open(fname, std::ofstream::out | std::ofstream::app);
		} else {
			f.open(fname);
		};

		// Make sure we found it
		if (!f.is_open()) {
			std::cerr << ">>> Error: Moment::write_to_file <<< could not write to file: " << fname << std::endl;
			exit(EXIT_FAILURE);
		};

		// Go through all time
        f << _val_averaged.at(MCType::AWAKE) << " " << _val_averaged.at(MCType::ASLEEP) << "\n";

		// Close
		f.close();
	};

};


