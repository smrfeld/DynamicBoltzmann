#include "../include/bmla_bits/ixn_param.hpp"

// Other headers
#include "../include/bmla_bits/general.hpp"
#include "../include/bmla_bits/species.hpp"
#include "../include/bmla_bits/moment.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

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
		double _lambda_s, _lambda_sp1;
		double _y_s, _y_sp1;

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

		void update_calculate_and_store(double dopt, bool l2_mode=false, double l2_lambda=0.0, double l2_center=0.0);
		void update_committ_stored(bool nesterov_mode=true);

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
		double lambda_0 = 0.0;
		_lambda_s = (1.0 + sqrt(1.0 + 4.0 * pow(lambda_0,2))) / 2.0;
		_lambda_sp1 = (1.0 + sqrt(1.0 + 4.0 * pow(_lambda_s,2))) / 2.0;
		_y_s = _val;
		_y_sp1 = _val;

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
		// ...
	};
	void IxnParam::Impl::_copy(const Impl& other) {
		_val = other._val;
		_update = other._update;
		_lambda_s = other._lambda_s;
		_lambda_sp1 = other._lambda_sp1;
		_y_s = other._y_s;
		_y_sp1 = other._y_sp1;

		_moment = other._moment;

		_is_val_fixed = other._is_val_fixed;
	};
	void IxnParam::Impl::_move(Impl& other) {
		_val = other._val;
		_update = other._update;
		_lambda_s = other._lambda_s;
		_lambda_sp1 = other._lambda_sp1;
		_y_s = other._y_s;
		_y_sp1 = other._y_sp1;

		_moment = std::move(other._moment);

		_is_val_fixed = other._is_val_fixed;

		// Reset the other
		other._val = 0.0;
		other._update = 0.0;
		other._lambda_s = 0.0;
		other._lambda_sp1 = 0.0;
		other._y_s = 0.0;
		other._y_sp1 = 0.0;
		other._is_val_fixed = false;
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
		return _val;
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

	void IxnParam::Impl::update_calculate_and_store(double dopt, bool l2_mode, double l2_lambda, double l2_center) {
		_update = - dopt * (_moment->get_moment(MomentType::ASLEEP) - _moment->get_moment(MomentType::AWAKE));
		if (l2_mode) {
			_update -= 2.0 * l2_lambda * (_val - l2_center);
		};
	};
	void IxnParam::Impl::update_committ_stored(bool nesterov_mode) {

		if (nesterov_mode) {

			// ysp1, lambda sp1
			_y_sp1 = _val + _update;
			// std::cout << "_y_sp1 = " << _y_sp1 << std::endl;
			_lambda_sp1 = (1.0 + sqrt(1.0 + 4.0 * pow(_lambda_s,2))) / 2.0;
			//std::cout << "_lambda_sp1 = " << _lambda_sp1 << std::endl;

			// Gamma
			double gamma_s = (1.0 - _lambda_s) / _lambda_sp1;
			// std::cout << "gamma_s = " << gamma_s << std::endl;

			// New val
			_val = (1.0 - gamma_s) * _y_sp1 + gamma_s * _y_s;

			// Advance y, lambda
			_y_s = _y_sp1;
			_lambda_s = _lambda_sp1;

		} else {

			// Just update
			_val += _update;

		};
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

	void IxnParam::update_calculate_and_store(double dopt, bool l2_mode, double l2_lambda, double l2_center) {
		_impl->update_calculate_and_store(dopt,l2_mode,l2_lambda,l2_center);
	};
	void IxnParam::update_committ_stored(bool nesterov_mode) {
		_impl->update_committ_stored(nesterov_mode);
	};

	/********************
	Write to file
	********************/

	void IxnParam::write_to_file(std::string fname, bool append) const {
		_impl->write_to_file(fname, append);
	};

};





