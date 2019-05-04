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
	Moment
	****************************************/
	
	/********************
	Constructor
	********************/

	Moment::Moment(std::string name, IxnParamType type) {
		_name = name;
		_type = type;

		// Fixed awake moment
		_is_awake_moment_fixed = false;
        
        // Offset
        _val_diff_offset = 0.0;
        
        // Vals
        _val_averaged[MCType::AWAKE] = 0.0;
        _val_averaged[MCType::ASLEEP] = 0.0;

        // Matrix of vals
        _weight_matrix[MCType::AWAKE] = nullptr;
        _weight_matrix[MCType::ASLEEP] = nullptr;
        _weight_matrix_awake_minus_asleep = nullptr;
        _bias_vec[MCType::AWAKE] = nullptr;
        _bias_vec[MCType::ASLEEP] = nullptr;
        _bias_vec_awake_minus_asleep = nullptr;
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

	void Moment::_clean_up() {
        if (_weight_matrix[MCType::AWAKE]) {
            delete _weight_matrix[MCType::AWAKE];
        };
        _weight_matrix[MCType::AWAKE] = nullptr;
        if (_weight_matrix[MCType::ASLEEP]) {
            delete _weight_matrix[MCType::ASLEEP];
        };
        _weight_matrix[MCType::ASLEEP] = nullptr;
        if (_weight_matrix_awake_minus_asleep) {
            delete _weight_matrix_awake_minus_asleep;
        };
        _weight_matrix_awake_minus_asleep = nullptr;
        
        if (_bias_vec[MCType::AWAKE]) {
            delete _bias_vec[MCType::AWAKE];
        };
        _bias_vec[MCType::AWAKE] = nullptr;
        if (_bias_vec[MCType::ASLEEP]) {
            delete _bias_vec[MCType::ASLEEP];
        };
        _bias_vec[MCType::ASLEEP] = nullptr;
        if (_bias_vec_awake_minus_asleep) {
            delete _bias_vec_awake_minus_asleep;
        };
        _bias_vec_awake_minus_asleep = nullptr;
    };
	void Moment::_move(Moment &other) {
		_name = other._name;
		_type = other._type;

		// averaged
		_val_averaged = other._val_averaged;
        _weight_matrix = other._weight_matrix;
        _weight_matrix_awake_minus_asleep = other._weight_matrix_awake_minus_asleep;
        _bias_vec = other._bias_vec;
        _bias_vec_awake_minus_asleep = other._bias_vec_awake_minus_asleep;
        
        _val_diff_offset = other._val_diff_offset;
        
		_is_awake_moment_fixed = other._is_awake_moment_fixed;

		// Reset the other
		other._name = "";

		other._val_averaged[MCType::AWAKE] = 0.0;
        other._val_averaged[MCType::ASLEEP] = 0.0;

        other._weight_matrix[MCType::AWAKE] = nullptr;
        other._weight_matrix[MCType::ASLEEP] = nullptr;
        other._weight_matrix_awake_minus_asleep = nullptr;
        
        other._bias_vec[MCType::AWAKE] = nullptr;
        other._bias_vec[MCType::ASLEEP] = nullptr;
        other._bias_vec_awake_minus_asleep = nullptr;
        
        other._val_diff_offset = 0.0;
        
		other._is_awake_moment_fixed = false;
	};
	void Moment::_copy(const Moment& other) {
		_name = other._name;
		_type = other._type;
        
        // _no_markov_chains = other._no_markov_chains;
        
		// reaped
        // _vals_reaped = other._vals_reaped;

		// averaged
		_val_averaged = other._val_averaged;
        
        if (other._weight_matrix.at(MCType::AWAKE)) {
            _weight_matrix[MCType::AWAKE] = new arma::sp_mat(*other._weight_matrix.at(MCType::AWAKE));
        } else {
            _weight_matrix[MCType::AWAKE] = nullptr;
        };
        if (other._weight_matrix.at(MCType::ASLEEP)) {
            _weight_matrix[MCType::ASLEEP] = new arma::sp_mat(*other._weight_matrix.at(MCType::ASLEEP));
        } else {
            _weight_matrix[MCType::ASLEEP] = nullptr;
        };
        if (other._weight_matrix_awake_minus_asleep) {
            _weight_matrix_awake_minus_asleep = new arma::sp_mat(*other._weight_matrix_awake_minus_asleep);
        } else {
            _weight_matrix_awake_minus_asleep = nullptr;
        };

        if (other._bias_vec.at(MCType::AWAKE)) {
            _bias_vec[MCType::AWAKE] = new arma::vec(*other._bias_vec.at(MCType::AWAKE));
        } else {
            _bias_vec[MCType::AWAKE] = nullptr;
        };
        if (other._bias_vec.at(MCType::ASLEEP)) {
            _bias_vec[MCType::ASLEEP] = new arma::vec(*other._bias_vec.at(MCType::ASLEEP));
        } else {
            _bias_vec[MCType::ASLEEP] = nullptr;
        };
        if (other._bias_vec_awake_minus_asleep) {
            _bias_vec_awake_minus_asleep = new arma::vec(*other._bias_vec_awake_minus_asleep);
        } else {
            _bias_vec_awake_minus_asleep = nullptr;
        };
        
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
    void Moment::increment_moment(MCType type, double val) {
        _val_averaged[type] += val;
    };
    void Moment::set_moment(MCType type, double val) {
        _val_averaged[type] = val;
    };
    void Moment::reset_moment(MCType type) {
        _val_averaged[type] = 0.0;
    };

    // Augment moment difference by some value
    void Moment::set_moment_offset(double val) {
        _val_diff_offset = val;
    };

    // Get moment difference
    double Moment::get_moment_diff_awake_minus_asleep_plus_offset() const {
        return _val_averaged.at(MCType::AWAKE) - _val_averaged.at(MCType::ASLEEP) + _val_diff_offset;
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

        f << _val_averaged.at(MCType::AWAKE) << " " << _val_averaged.at(MCType::ASLEEP) << "\n";

		// Close
		f.close();
	};

    void Moment::write_weight_matrix_to_file(std::string fname) const {
        std::ofstream f;
        
        // Open
        f.open(fname);
        
        // Make sure we found it
        if (!f.is_open()) {
            std::cerr << ">>> Error: Moment::write_weight_matrix_to_file <<< could not write to file: " << fname << std::endl;
            exit(EXIT_FAILURE);
        };
        
        for (auto i=0; i<_weight_matrix.at(MCType::AWAKE)->n_rows; i++) {
            for (auto j=0; j<_weight_matrix.at(MCType::AWAKE)->n_cols; j++) {
                if ((*_weight_matrix_awake_minus_asleep)(i,j) != 0) {
                    f << i << " " << j << " " << (*_weight_matrix_awake_minus_asleep)(i,j) << "\n";
                };
            };
        };
        
        // Close
        f.close();
    };
};


