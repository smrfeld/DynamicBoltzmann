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
            _weight_matrix[MCType::AWAKE] = new arma::mat(*other._weight_matrix.at(MCType::AWAKE));
        } else {
            _weight_matrix[MCType::AWAKE] = nullptr;
        };
        if (other._weight_matrix.at(MCType::ASLEEP)) {
            _weight_matrix[MCType::ASLEEP] = new arma::mat(*other._weight_matrix.at(MCType::ASLEEP));
        } else {
            _weight_matrix[MCType::ASLEEP] = nullptr;
        };
        if (other._weight_matrix_awake_minus_asleep) {
            _weight_matrix_awake_minus_asleep = new arma::mat(*other._weight_matrix_awake_minus_asleep);
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
    
    // ***************
    // MARK: - Set W matrix / bias vec
    // ***************
    
    void Moment::set_weight_matrix_dims(int dim_upper_layer, int dim_lower_layer) {
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
        
        if (_type != IxnParamType::W && _type != IxnParamType::X) {
            std::cerr << ">>> Moment::set_weight_matrix_dims <<< only W and X ixn params can have a weight matrix" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        _weight_matrix[MCType::AWAKE] = new arma::mat(dim_upper_layer,dim_lower_layer,arma::fill::zeros);
        _weight_matrix[MCType::ASLEEP] = new arma::mat(dim_upper_layer,dim_lower_layer,arma::fill::zeros);
        _weight_matrix_awake_minus_asleep = new arma::mat(dim_upper_layer,dim_lower_layer,arma::fill::zeros);
    };
    void Moment::reset_weight_matrix(MCType type) {
        _weight_matrix[type]->fill(arma::fill::zeros);
    };
    void Moment::increment_weight_matrix(MCType type, arma::mat increment) {
        *_weight_matrix[type] += increment;
    };
    void Moment::set_moment_to_weight_matrix_sum(MCType type) {
        _val_averaged[type] = arma::accu(*_weight_matrix[type]);
    };
    const arma::mat& Moment::get_weight_matrix(MCType type) const {
        return *_weight_matrix.at(type);
    };

    void Moment::calculate_weight_matrix_awake_minus_asleep() {
        *_weight_matrix_awake_minus_asleep = *_weight_matrix.at(MCType::AWAKE) - *_weight_matrix.at(MCType::ASLEEP);
    };
    const arma::mat& Moment::get_weight_matrix_awake_minus_asleep() const {
        return *_weight_matrix_awake_minus_asleep;
    };
    
    void Moment::set_bias_vec_dims(int dims) {
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
        
        if (_type != IxnParamType::B && _type != IxnParamType::H) {
            std::cerr << ">>> Moment::set_bias_vec_dims <<< only H and B ixn params can have a bias vec" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        _bias_vec[MCType::AWAKE] = new arma::vec(dims,arma::fill::zeros);
        _bias_vec[MCType::ASLEEP] = new arma::vec(dims,arma::fill::zeros);
        _bias_vec_awake_minus_asleep = new arma::vec(dims,arma::fill::zeros);
    };
    void Moment::reset_bias_vec(MCType type) {
        _bias_vec[type]->fill(arma::fill::zeros);
    };
    void Moment::increment_bias_vec(MCType type, arma::vec increment) {
        *_bias_vec[type] += increment;
    };
    void Moment::set_moment_to_bias_vec_sum(MCType type) {
        _val_averaged[type] = arma::accu(*_bias_vec.at(type));
    };
    const arma::vec& Moment::get_bias_vec(MCType type) const {
        return *_bias_vec.at(type);
    };

    void Moment::calculate_bias_vec_awake_minus_asleep() {
        *_bias_vec_awake_minus_asleep = *_bias_vec.at(MCType::AWAKE) - *_bias_vec.at(MCType::ASLEEP);
    };
    const arma::vec& Moment::get_bias_vec_awake_minus_asleep() const {
        return *_bias_vec_awake_minus_asleep;
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


