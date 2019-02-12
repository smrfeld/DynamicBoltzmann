#include "../include/dblz_bits/ixn_param.hpp"

// Other headers
#include "../include/dblz_bits/general.hpp"
#include "../include/dblz_bits/species.hpp"
#include "../include/dblz_bits/moment.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

/************************************
* Namespace for bmla
************************************/

namespace dblz {

	/****************************************
	IxnFunc
	****************************************/

    // ***************
    // MARK: - Constructor
    // ***************
    
	IxnParam::IxnParam(std::string name, IxnParamType type, double init_guess) {

		_val = init_guess;
		_update = 0.0;
        
		// nesterov
		// double lambda_0 = 0.0;
		// _lambda_s = (1.0 + sqrt(1.0 + 4.0 * pow(lambda_0,2))) / 2.0;
		// _lambda_sp1 = (1.0 + sqrt(1.0 + 4.0 * pow(_lambda_s,2))) / 2.0;
		_nesterov_y_s = nullptr;
		_nesterov_y_sp1 = nullptr;

		// adam
		_adam_v = nullptr;
		_adam_m = nullptr;

		// Moment
		_moment = std::make_unique<Moment>(name, type);

		_is_val_fixed = false;
	};
	IxnParam::IxnParam(const IxnParam& other) {
		_copy(other);
	};
	IxnParam::IxnParam(IxnParam&& other) {
		_move(other);
	};
    IxnParam& IxnParam::operator=(const IxnParam& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    IxnParam& IxnParam::operator=(IxnParam&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };
	IxnParam::~IxnParam()
	{
		_clean_up();
	};
	void IxnParam::_clean_up() {
		if (_nesterov_y_s) {
			delete _nesterov_y_s;
		};
        _nesterov_y_s = nullptr;
        
        if (_nesterov_y_sp1) {
			delete _nesterov_y_sp1;
		};
        _nesterov_y_sp1 = nullptr;

        if (_adam_m) {
			delete _adam_m;
		};
        _adam_m = nullptr;

        if (_adam_v) {
			delete _adam_v;
		};
        _adam_v = nullptr;
    };
	void IxnParam::_copy(const IxnParam& other) {
		_val = other._val;
		_update = other._update;
        if (other._nesterov_y_s) {
            _nesterov_y_s = new double(*other._nesterov_y_s);
        } else {
            _nesterov_y_s = nullptr;
        };
        if (other._nesterov_y_sp1) {
            _nesterov_y_sp1 = new double(*other._nesterov_y_sp1);
        } else {
            _nesterov_y_sp1 = nullptr;
        };
        if (other._adam_m) {
            _adam_m = new double(*other._adam_m);
        } else {
            _adam_m = nullptr;
        };
        if (other._adam_v) {
            _adam_v = new double(*other._adam_v);
        } else {
            _adam_v = nullptr;
        };

        _moment = std::make_shared<Moment>(*other._moment);

		_is_val_fixed = other._is_val_fixed;
	};
	void IxnParam::_move(IxnParam& other) {
		_val = other._val;
		_update = other._update;
        _nesterov_y_s = other._nesterov_y_s;
        _nesterov_y_sp1 = other._nesterov_y_sp1;
        _adam_m = other._adam_m;
        _adam_v = other._adam_v;
        
		_moment = std::move(other._moment);

		_is_val_fixed = other._is_val_fixed;

		// Reset the other
		other._val = 0.0;
		other._update = 0.0;
		other._nesterov_y_s = nullptr;
		other._nesterov_y_sp1 = nullptr;
		other._is_val_fixed = false;
		other._adam_m = nullptr;
		other._adam_v = nullptr;
	};

	/********************
	Name, type
	********************/

	std::string IxnParam::get_name() const {
		return _moment->get_name();
	};

	IxnParamType IxnParam::get_type() const {
		return _moment->get_type();
	};

	/********************
	Value
	********************/

	double IxnParam::get_val() const {
		return  _val;
	};
	void IxnParam::set_val(double val) {
		_val = val;
	};

	/********************
	Fixed value
	********************/

	void IxnParam::set_fix_value(bool fixed) {
		_is_val_fixed = fixed;
	};
	bool IxnParam::get_is_val_fixed() const {
		return _is_val_fixed;
	};

	/********************
	Moment
	********************/

	std::shared_ptr<Moment> IxnParam::get_moment() const {
		return _moment;
	};

	/********************
	Update
	********************/

	void IxnParam::update_calculate_and_store_l2(double l2_lambda, double l2_center) {
		_update = _moment->get_moment_diff_awake_minus_asleep() - 2.0 * l2_lambda * (_val - l2_center);
	};
    void IxnParam::update_calculate_and_store() {
        _update = _moment->get_moment_diff_awake_minus_asleep();
    };
    
    void IxnParam::update_committ_stored_sgd(double dopt) {
		// Just update
		_val += dopt*_update;
	};
	void IxnParam::update_committ_stored_nesterov(double dopt, double nesterov_acc) {

		// ysp1, lambda sp1
		if (!_nesterov_y_sp1) {
			_nesterov_y_sp1 = new double(_val);
		};
		*_nesterov_y_sp1 = _val - dopt * _update;
		// std::cout << "_nesterov_y_sp1 = " << _nesterov_y_sp1 << std::endl;
		// _lambda_sp1 = (1.0 + sqrt(1.0 + 4.0 * pow(_lambda_s,2))) / 2.0;
		//std::cout << "_lambda_sp1 = " << _lambda_sp1 << std::endl;

		// Gamma
		// double gamma_s = (1.0 - _lambda_s) / _lambda_sp1;
		// std::cout << "gamma_s = " << gamma_s << std::endl;

		// Y_S
		if (!_nesterov_y_s) {
			_nesterov_y_s = new double(_val);
		};

		// New val
		// _val = (1.0 - gamma_s) * _nesterov_y_sp1 + gamma_s * _nesterov_y_s;
		_val = (1.0 + nesterov_acc) * (*_nesterov_y_sp1) - nesterov_acc * (*_nesterov_y_s);

		// Advance y, lambda
		*_nesterov_y_s = *_nesterov_y_sp1;
		// _lambda_s = _lambda_sp1;
	};
	void IxnParam::update_committ_stored_adam(double dopt, int opt_step, double beta_1, double beta_2, double eps) {
		if (!_adam_m) {
			_adam_m = new double(0.0);
		};
		if (!_adam_v) {
			_adam_v = new double(0.0);
		};

		int opt_step_use = std::max(opt_step,1);

		// first moment
		*_adam_m = beta_1*(*_adam_m) + (1.0 - beta_1)*_update;

		// second
		*_adam_v = beta_2*(*_adam_v) + (1.0-beta_2)*pow(_update,2);

		// corrections
		double mhat = (*_adam_m) / (1.0 - pow(beta_1,opt_step_use));
		double vhat = (*_adam_v) / (1.0 - pow(beta_2,opt_step_use));

		// update
		_val += dopt * mhat / (sqrt(vhat) + eps);
	};

	/********************
	Write to file
	********************/

	void IxnParam::write_to_file(std::string fname, bool append) const {
		std::ofstream f;

		// Open
		if (append) {
			f.open(fname, std::ofstream::out | std::ofstream::app);
		} else {
			f.open(fname);
		};

		// Make sure we found it
		if (!f.is_open()) {
			std::cerr << ">>> Error: IxnParam::write_to_file <<< could not write to file: " << fname << std::endl;
			exit(EXIT_FAILURE);
		};

		// Go through all time
		f << _val << "\n";

		// Close
		f.close();
	};
};





