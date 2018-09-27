#include "../include/bmla_bits/ixn_param.hpp"

// Other headers
#include "../include/bmla_bits/general.hpp"
#include "../include/bmla_bits/species.hpp"
#include "../include/bmla_bits/moment.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

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

		/********************
		Moment
		********************/

		std::shared_ptr<Moment> get_moment() const;

		/********************
		Update
		********************/

		void update_calculate_and_store(double dopt);
		void update_committ_stored();

		/********************
		Write to file
		********************/

		void write_to_file(std::string fname) const;

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

		// Moment
		_moment = std::make_unique<Moment>(name, type);
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

		_moment = other._moment;
	};
	void IxnParam::Impl::_move(Impl& other) {
		_val = other._val;
		_update = other._update;

		_moment = std::move(other._moment);

		// Reset the other
		other._val = 0.0;
		other._update = 0.0;
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

	/********************
	Moment
	********************/

	std::shared_ptr<Moment> IxnParam::Impl::get_moment() const {
		return _moment;
	};

	/********************
	Update
	********************/

	void IxnParam::Impl::update_calculate_and_store(double dopt) {
		_update = - dopt * (_moment->get_moment(MomentType::ASLEEP) - _moment->get_moment(MomentType::AWAKE));
	};
	void IxnParam::Impl::update_committ_stored() {
		_val += _update;
	};

	/********************
	Write to file
	********************/

	void IxnParam::Impl::write_to_file(std::string fname) const {
		std::ofstream f;

		// Open
		f.open(fname);

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
	Value
	********************/

	double IxnParam::get_val() const {
		return _impl->get_val();
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

	void IxnParam::update_calculate_and_store(double dopt) {
		_impl->update_calculate_and_store(dopt);
	};
	void IxnParam::update_committ_stored() {
		_impl->update_committ_stored();
	};

	/********************
	Write to file
	********************/

	void IxnParam::write_to_file(std::string fname) const {
		_impl->write_to_file(fname);
	};

};





