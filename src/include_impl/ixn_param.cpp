#include "../../include/dynamicboltz_bits/ixn_param.hpp"

// Other headers
#include "../../include/dynamicboltz_bits/general.hpp"
#include "../../include/dynamicboltz_bits/species.hpp"
#include "../../include/dynamicboltz_bits/moment.hpp"

#include <iostream>
#include <fstream>

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	Ixn Param Traj Implementation
	****************************************/

	class IxnParam::Impl {

	private:

		// Values
		// Timepoints = timesteps + 1
		int _no_timepoints;
		int _no_timesteps;
		double *_vals;
		double _init_cond;

		// Diff eq RHS
		std::shared_ptr<DiffEqRHS> _diff_eq;

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

		Impl(std::string name, IxnParamType type, double init_cond); 
		Impl(const Impl& other);
		Impl(Impl&& other);
		Impl& operator=(const Impl& other);
		Impl& operator=(Impl&& other);
		~Impl();

		/********************
		Validate setup
		********************/

		void check_setup() const;

		/********************
		Timepoints
		********************/

		int get_no_timesteps() const;
		void set_no_timesteps(int no_timesteps);

		/********************
		Init cond
		********************/

		double get_init_cond() const;
		void set_init_cond(double init_cond);

		// void set_fixed_awake_moment(std::vector<double> vals);

		/********************
		Name, type
		********************/

		std::string get_name() const;

		IxnParamType get_type() const;

		/********************
		Value
		********************/

		double get_val_at_timepoint(int timepoint) const;

		/********************
		Diff eq
		********************/

		void set_diff_eq(std::shared_ptr<DiffEqRHS> diff_eq);

		void solve_diff_eq_at_timepoint_to_plus_one(int timepoint, double dt);

		/********************
		Moment
		********************/

		std::shared_ptr<Moment> get_moment() const;

	};








































	/****************************************
	IxnFunc - IMPLEMENTATION DEFINITIONS
	****************************************/

	/********************
	Constructor
	********************/

	IxnParam::Impl::Impl(std::string name, IxnParamType type, double init_cond) {

		_init_cond = init_cond;

		_no_timesteps = 0;
		_no_timepoints = 1;

		_vals = new double[_no_timepoints];
		std::fill_n(_vals,_no_timepoints,0.0);
		_vals[0] = _init_cond;

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
		safeDelArr(_vals);
	};
	void IxnParam::Impl::_copy(const Impl& other) {
		_no_timesteps = other._no_timesteps;
		_no_timepoints = other._no_timepoints;
		_vals = new double[_no_timepoints];
		std::copy(other._vals,other._vals+_no_timepoints,_vals);
		_init_cond = other._init_cond;

		_diff_eq = other._diff_eq;
		_moment = other._moment;
	};
	void IxnParam::Impl::_move(Impl& other) {
		_no_timesteps = other._no_timesteps;
		_no_timepoints = other._no_timepoints;
		_vals = new double[_no_timepoints];
		std::copy(other._vals,other._vals+_no_timepoints,_vals);
		_init_cond = other._init_cond;

		_diff_eq = std::move(other._diff_eq);
		_moment = std::move(other._moment);

		// Reset the other
		other._no_timesteps = 0;
		other._no_timepoints = 0;
		safeDelArr(other._vals);
		other._init_cond = 0.0;
	};

	/********************
	Validate setup
	********************/

	void IxnParam::Impl::check_setup() const {
		// ...
	};

	/********************
	Timesteps
	********************/

	int IxnParam::Impl::get_no_timesteps() const {
		return _no_timesteps;
	};
	void IxnParam::Impl::set_no_timesteps(int no_timesteps) {
		_no_timesteps = no_timesteps;
		_no_timepoints = _no_timesteps + 1;
		// Vals
		safeDelArr(_vals);
		_vals = new double[_no_timepoints];
		std::fill_n(_vals,_no_timepoints,0.0);
		_vals[0] = _init_cond;
	};


	/********************
	Init cond
	********************/

	double IxnParam::Impl::get_init_cond() const {
		return _init_cond;
	};
	void IxnParam::Impl::set_init_cond(double init_cond) {
		_init_cond = init_cond;
		_vals[0] = _init_cond;
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

	double IxnParam::Impl::get_val_at_timepoint(int timepoint) const {
		if (timepoint >= _no_timepoints) {
			std::cerr << ">>> Error: IxnParam::Impl::get_val_at_timepoint <<< " << timepoint << " is out of bounds: " << _no_timepoints << std::endl;
		};
		return _vals[timepoint];
	};

	/********************
	Diff eq
	********************/

	void IxnParam::Impl::set_diff_eq(std::shared_ptr<DiffEqRHS> diff_eq) {
		_diff_eq = diff_eq;
	};

	void IxnParam::Impl::solve_diff_eq_at_timepoint_to_plus_one(int timepoint, double dt) {
		// ...
	};

	/********************
	Moment
	********************/

	std::shared_ptr<Moment> IxnParam::Impl::get_moment() const {
		return _moment;
	};





































	/****************************************
	IxnFuncSpec - Impl forwards
	****************************************/

	/********************
	Constructor
	********************/

	IxnParam::IxnParam(std::string name, IxnParamType type, double init_cond) : _impl(new Impl(name,type,init_cond)) {};
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
	Validate setup
	********************/

	void IxnParam::check_setup() const {
		_impl->check_setup();
	};

	/********************
	Timesteps
	********************/

	int IxnParam::get_no_timesteps() const {
		return _impl->get_no_timesteps();
	};
	void IxnParam::set_no_timesteps(int no_timesteps) {
		_impl->set_no_timesteps(no_timesteps);
	};

	/********************
	Init cond
	********************/

	double IxnParam::get_init_cond() const {
		return _impl->get_init_cond();
	};
	void IxnParam::set_init_cond(double init_cond) {
		_impl->set_init_cond(init_cond);
	};

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

	double IxnParam::get_val_at_timepoint(int timepoint) const {
		return _impl->get_val_at_timepoint(timepoint);
	};

	/********************
	Diff eq
	********************/

	void IxnParam::set_diff_eq(std::shared_ptr<DiffEqRHS> diff_eq) {
		_impl->set_diff_eq(diff_eq);
	};

	void IxnParam::solve_diff_eq_at_timepoint_to_plus_one(int timepoint, double dt) {
		_impl->solve_diff_eq_at_timepoint_to_plus_one(timepoint,dt);
	};

	/********************
	Moment
	********************/

	std::shared_ptr<Moment> IxnParam::get_moment() const {
		return _impl->get_moment();
	};
};





