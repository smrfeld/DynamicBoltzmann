#include "../include/bmla_bits/ixn_param.hpp"

// Other headers
#include "../include/bmla_bits/general.hpp"
#include "../include/bmla_bits/species.hpp"
#include "../include/bmla_bits/moment.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

/************************************
* Namespace for bmla
************************************/

namespace bmla {

	/****************************************
	Ixn Param Traj Implementation
	****************************************/

	class IxnParam::Impl {

	private:

		// Values
		double _val;
		double _update;
        
		// nesterov
		// double _lambda_s, _lambda_sp1;
		double *_nesterov_y_s, *_nesterov_y_sp1;

		// adam
		double *_adam_m;
		double *_adam_v;

		// Fixed value
		bool _is_val_fixed;

		// Moment
		std::shared_ptr<Moment> _moment;

		// Copy, clean up
		void _clean_up();
		void _copy(const Impl& other);
		void _move(Impl &other);

	public:

		/********************
		Constructor
		********************/

		Impl(std::string name, IxnParamType type, double init_guess); 
		Impl(const Impl& other);
		Impl(Impl&& other);
		Impl& operator=(const Impl& other);
		Impl& operator=(Impl&& other);
		~Impl();

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

		void update_calculate_and_store(bool l2_mode=false, double l2_lambda=0.0, double l2_center=0.0);
		void update_committ_stored_sgd(double dopt);
		void update_committ_stored_nesterov(double dopt, double nesterov_acc);
		void update_committ_stored_adam(double dopt, int opt_step, double beta_1, double beta_2, double eps);

		/********************
		Write to file
		********************/

		void write_to_file(std::string fname, bool append=false) const;

	};








































	/****************************************
	IxnFunc - IMPLEMENTATION DEFINITIONS
	****************************************/

	/********************
	Constructor
	********************/

	IxnParam::Impl::Impl(std::string name, IxnParamType type, double init_guess) {

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
	IxnParam::Impl::Impl(const Impl& other) {
		_copy(other);
	};
	IxnParam::Impl::Impl(Impl&& other) {
		_move(other);
	};
    IxnParam::Impl& IxnParam::Impl::operator=(const Impl& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    IxnParam::Impl& IxnParam::Impl::operator=(Impl&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };
	IxnParam::Impl::~Impl()
	{
		_clean_up();
	};
	void IxnParam::Impl::_clean_up() {
		if (_nesterov_y_s) {
			delete _nesterov_y_s;
			_nesterov_y_s = nullptr;
		};
		if (_nesterov_y_sp1) {
			delete _nesterov_y_sp1;
			_nesterov_y_sp1 = nullptr;
		};
		if (_adam_m) {
			delete _adam_m;
			_adam_m = nullptr;
		};
		if (_adam_v) {
			delete _adam_v;
			_adam_v = nullptr;
		};	
	};
	void IxnParam::Impl::_copy(const Impl& other) {
		_val = other._val;
		_update = other._update;
		// _lambda_s = other._lambda_s;
		// _lambda_sp1 = other._lambda_sp1;
		_nesterov_y_s = new double(*other._nesterov_y_s);
		_nesterov_y_sp1 = new double(*other._nesterov_y_sp1);
		_adam_m = new double(*other._adam_m);
		_adam_v = new double(*other._adam_v);

		_moment = other._moment;

		_is_val_fixed = other._is_val_fixed;
	};
	void IxnParam::Impl::_move(Impl& other) {
		_val = other._val;
		_update = other._update;
		// _lambda_s = other._lambda_s;
		// _lambda_sp1 = other._lambda_sp1;
		_nesterov_y_s = other._nesterov_y_s;
		_nesterov_y_sp1 = other._nesterov_y_sp1;
		_adam_m = other._adam_m;
		_adam_v = other._adam_v;

		_moment = std::move(other._moment);

		_is_val_fixed = other._is_val_fixed;

		// Reset the other
		other._val = 0.0;
		other._update = 0.0;
		// other._lambda_s = 0.0;
		// other._lambda_sp1 = 0.0;
		other._nesterov_y_s = nullptr;
		other._nesterov_y_sp1 = nullptr;
		other._is_val_fixed = false;
		other._adam_m = nullptr;
		other._adam_v = nullptr;
	};

	/********************
	Name, type
	********************/

	std::string IxnParam::Impl::get_name() const {
		return _moment->get_name();
	};

	IxnParamType IxnParam::Impl::get_type() const {
		return _moment->get_type();
	};

	/********************
	Value
	********************/

	double IxnParam::Impl::get_val() const {
		return  _val;
	};
	void IxnParam::Impl::set_val(double val) {
		_val = val;
	};

	/********************
	Fixed value
	********************/

	void IxnParam::Impl::set_fix_value(bool fixed) {
		_is_val_fixed = fixed;
	};
	bool IxnParam::Impl::get_is_val_fixed() const {
		return _is_val_fixed;
	};

	/********************
	Moment
	********************/

	std::shared_ptr<Moment> IxnParam::Impl::get_moment() const {
		return _moment;
	};

	/********************
	Update
	********************/

	void IxnParam::Impl::update_calculate_and_store(bool l2_mode, double l2_lambda, double l2_center) {
		_update = _moment->get_moment(MomentType::ASLEEP) - _moment->get_moment(MomentType::AWAKE);
		if (l2_mode) {
			_update += 2.0 * l2_lambda * (_val - l2_center);
		};
	};
	void IxnParam::Impl::update_committ_stored_sgd(double dopt) {
		// Just update
		_val -= dopt*_update;
	};
	void IxnParam::Impl::update_committ_stored_nesterov(double dopt, double nesterov_acc) {

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
	void IxnParam::Impl::update_committ_stored_adam(double dopt, int opt_step, double beta_1, double beta_2, double eps) {
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
		_val -= dopt * mhat / (sqrt(vhat) + eps);
	};

	/********************
	Write to file
	********************/

	void IxnParam::Impl::write_to_file(std::string fname, bool append) const {
		std::ofstream f;

		// Open
		if (append) {
			f.open(fname, std::ofstream::out | std::ofstream::app);
		} else {
			f.open(fname);
		};

		// Make sure we found it
		if (!f.is_open()) {
			std::cerr << ">>> Error: IxnParam::Impl::write_to_file <<< could not write to file: " << fname << std::endl;
			exit(EXIT_FAILURE);
		};

		// Go through all time
		f << _val << "\n";

		// Close
		f.close();
	};
































	/****************************************
	IxnFuncSpec - Impl forwards
	****************************************/

	/********************
	Constructor
	********************/

	IxnParam::IxnParam(std::string name, IxnParamType type, double init_guess) : _impl(new Impl(name,type,init_guess)) {};
	IxnParam::IxnParam(const IxnParam& other) : _impl(new Impl(*other._impl)) {};
	IxnParam::IxnParam(IxnParam&& other) : _impl(std::move(other._impl)) {};
	IxnParam& IxnParam::operator=(const IxnParam &other) {
        _impl.reset( new Impl( *other._impl ) );
        return *this; 
	};
	IxnParam& IxnParam::operator=(IxnParam &&other) {
        _impl = std::move(other._impl);
        return *this; 
	};
	IxnParam::~IxnParam() = default;

	/********************
	Name, type
	********************/

	std::string IxnParam::get_name() const {
		return _impl->get_name();
	};

	IxnParamType IxnParam::get_type() const {
		return _impl->get_type();
	};

	/********************
	Fixed value
	********************/

	void IxnParam::set_fix_value(bool fixed) {
		_impl->set_fix_value(fixed);
	};
	bool IxnParam::get_is_val_fixed() const {
		return _impl->get_is_val_fixed();
	};

	/********************
	Value
	********************/

	double IxnParam::get_val() const {
		return _impl->get_val();
	};
	void IxnParam::set_val(double val) {
		_impl->set_val(val);
	};

	/********************
	Moment
	********************/

	std::shared_ptr<Moment> IxnParam::get_moment() const {
		return _impl->get_moment();
	};

	/********************
	Update
	********************/

	void IxnParam::update_calculate_and_store(bool l2_mode, double l2_lambda, double l2_center) {
		_impl->update_calculate_and_store(l2_mode,l2_lambda,l2_center);
	};
	void IxnParam::update_committ_stored_sgd(double dopt) {
		_impl->update_committ_stored_sgd(dopt);
	};
	void IxnParam::update_committ_stored_nesterov(double dopt, double nesterov_acc) {
		_impl->update_committ_stored_nesterov(dopt,nesterov_acc);
	};
	void IxnParam::update_committ_stored_adam(double dopt, int opt_step, double beta_1, double beta_2, double eps) {
		_impl->update_committ_stored_adam(dopt,opt_step,beta_1,beta_2,eps);
	};

	/********************
	Write to file
	********************/

	void IxnParam::write_to_file(std::string fname, bool append) const {
		_impl->write_to_file(fname, append);
	};

};





