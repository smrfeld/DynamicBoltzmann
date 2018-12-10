#include "../include/dynamicboltz_bits/ixn_param.hpp"

// Other headers
#include "../include/dynamicboltz_bits/general.hpp"
#include "../include/dynamicboltz_bits/species.hpp"
#include "../include/dynamicboltz_bits/moment.hpp"
#include "../include/dynamicboltz_bits/diff_eq_rhs.hpp"
#include "../include/dynamicboltz_bits/adjoint.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	Ixn Param Traj Implementation
	****************************************/

	class IxnParam::Impl {

	private:

		// Adjoint
		std::shared_ptr<Adjoint> _adjoint;

		// Values
		// Timepoints = timesteps + 1
		int _no_timepoints;
		int _no_timesteps;
		double *_vals;
		double _init_cond;

		// Fixed value at the init cond
		bool _is_val_fixed_to_init_cond;
		bool _are_vals_fixed;

		// Diff eq RHS
		std::shared_ptr<DiffEqRHS> _diff_eq;

		// Moment
		std::shared_ptr<Moment> _moment;

		// Diff eq rhs that this ixn param appears in
		// and the corresponding idx in that function (= 0, 1, 2, etc to denote argument order)
		std::vector<std::pair<std::shared_ptr<DiffEqRHS>,int>> _diff_eq_dependencies;

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
		Diff eq rhs that this ixn param appears in
		********************/

		void add_diff_eq_dependency(std::shared_ptr<DiffEqRHS> diff_eq, int idx);
		const std::vector<std::pair<std::shared_ptr<DiffEqRHS>,int>>& get_diff_eq_dependencies() const;

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

		/********************
		Fixed value to IC
		********************/

		void set_is_val_fixed_to_init_cond(bool fixed);
		bool get_is_val_fixed_to_init_cond() const;

		void set_are_vals_fixed(bool fixed);
		bool get_are_vals_fixed() const;

		/********************
		Name, type
		********************/

		std::string get_name() const;

		IxnParamType get_type() const;

		/********************
		Value
		********************/

		void set_val_at_timepoint(int timepoint, double val);
		double get_val_at_timepoint(int timepoint) const;

		/********************
		Diff eq
		********************/

		std::shared_ptr<DiffEqRHS> get_diff_eq_rhs() const;
		void set_diff_eq_rhs(std::shared_ptr<DiffEqRHS> diff_eq);

		void solve_diff_eq_at_timepoint_to_plus_one(int timepoint, double dt);

		/********************
		Moment
		********************/

		std::shared_ptr<Moment> get_moment() const;

		/********************
		Adjoint
		********************/

		void set_adjoint(std::shared_ptr<Adjoint> adjoint);
		std::shared_ptr<Adjoint> get_adjoint() const;

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

	IxnParam::Impl::Impl(std::string name, IxnParamType type, double init_cond) {

		_init_cond = init_cond;

		_no_timesteps = 0;
		_no_timepoints = 1;

		_vals = new double[_no_timepoints];
		std::fill_n(_vals,_no_timepoints,0.0);
		_vals[0] = _init_cond;

		// Moment
		_moment = std::make_unique<Moment>(name, type);

		// Adjoint
		_adjoint = nullptr;

		_is_val_fixed_to_init_cond = false;
		_are_vals_fixed = false;
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
		_diff_eq_dependencies = other._diff_eq_dependencies;
		_moment = other._moment;
		_adjoint = other._adjoint;
		_is_val_fixed_to_init_cond = other._is_val_fixed_to_init_cond;
		_are_vals_fixed = other._are_vals_fixed;
	};
	void IxnParam::Impl::_move(Impl& other) {
		_no_timesteps = other._no_timesteps;
		_no_timepoints = other._no_timepoints;
		_vals = other._vals;
		_init_cond = other._init_cond;

		_diff_eq = std::move(other._diff_eq);
		_diff_eq_dependencies = other._diff_eq_dependencies;
		_moment = std::move(other._moment);
		_is_val_fixed_to_init_cond = other._is_val_fixed_to_init_cond;
		_are_vals_fixed = other._are_vals_fixed;

		// Reset the other
		other._diff_eq_dependencies.clear();
		other._no_timesteps = 0;
		other._no_timepoints = 0;
		other._vals = nullptr;
		other._init_cond = 0.0;
		other._is_val_fixed_to_init_cond = false;
		other._are_vals_fixed = false;

		_adjoint = std::move(other._adjoint);
	};

	/********************
	Diff eq rhs that this ixn param appears in
	********************/

	void IxnParam::Impl::add_diff_eq_dependency(std::shared_ptr<DiffEqRHS> diff_eq, int idx) {
		_diff_eq_dependencies.push_back(std::make_pair(diff_eq,idx));
	};
	const std::vector<std::pair<std::shared_ptr<DiffEqRHS>,int>>& IxnParam::Impl::get_diff_eq_dependencies() const {
		return _diff_eq_dependencies;
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

		// If fixed, fill with IC
		set_is_val_fixed_to_init_cond(_is_val_fixed_to_init_cond);

		// If fixed, set val
		set_are_vals_fixed(_are_vals_fixed);

		// Set for moment
		_moment->set_no_timesteps(_no_timesteps);
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

		// update if fixed
		set_is_val_fixed_to_init_cond(_is_val_fixed_to_init_cond);
	};

	/********************
	Fixed value to IC
	********************/

	void IxnParam::Impl::set_is_val_fixed_to_init_cond(bool fixed) {
		_is_val_fixed_to_init_cond = fixed;
		if (_is_val_fixed_to_init_cond) {
			std::fill_n(_vals,_no_timepoints,_init_cond);
		};
	};
	bool IxnParam::Impl::get_is_val_fixed_to_init_cond() const {
		return _is_val_fixed_to_init_cond;
	};

	void IxnParam::Impl::set_are_vals_fixed(bool fixed) {
		_are_vals_fixed = fixed;
	};
	bool IxnParam::Impl::get_are_vals_fixed() const {
		return _are_vals_fixed;
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

	void IxnParam::Impl::set_val_at_timepoint(int timepoint, double val) {
		if (timepoint >= _no_timepoints) {
			std::cerr << ">>> Error: IxnParam::Impl::set_val_at_timepoint <<< " << timepoint << " is out of bounds: " << _no_timepoints << std::endl;
			exit(EXIT_FAILURE);
		};
		_vals[timepoint] = val;
	};
	double IxnParam::Impl::get_val_at_timepoint(int timepoint) const {
		if (timepoint >= _no_timepoints) {
			std::cerr << ">>> Error: IxnParam::Impl::get_val_at_timepoint <<< " << timepoint << " is out of bounds: " << _no_timepoints << std::endl;
			exit(EXIT_FAILURE);
		};
		return _vals[timepoint];
	};

	/********************
	Diff eq
	********************/

	std::shared_ptr<DiffEqRHS> IxnParam::Impl::get_diff_eq_rhs() const {
		return _diff_eq;
	};
	void IxnParam::Impl::set_diff_eq_rhs(std::shared_ptr<DiffEqRHS> diff_eq) {
		_diff_eq = diff_eq;
	};

	void IxnParam::Impl::solve_diff_eq_at_timepoint_to_plus_one(int timepoint, double dt) {
		if (_is_val_fixed_to_init_cond) { // Do nothing if val is fixed
			std::cerr << ">>> Error: IxnParam::Impl::solve_diff_eq_at_timepoint_to_plus_one <<< ixn param: " << get_name() << " is fixed!" << std::endl;
			exit(EXIT_FAILURE);
		};

		if (timepoint >= _no_timepoints) {
			std::cerr << ">>> Error: IxnParam::Impl::solve_diff_eq_at_timepoint_to_plus_one <<< " << timepoint << " is out of bounds: " << _no_timepoints << std::endl;
			exit(EXIT_FAILURE);
		};
		if (!_diff_eq) {
			std::cerr << ">>> Error: IxnParam::Impl::solve_diff_eq_at_timepoint_to_plus_one <<< No diff eq set" << std::endl;
			exit(EXIT_FAILURE);
		};

		// Check in domain
		if (!_diff_eq->check_val_is_in_domain_at_timepoint(timepoint)) {
			// std::cerr << ">>> Error: IxnParam::Impl::solve_diff_eq_at_timepoint_to_plus_one <<< outside of domain of ixn param: " << get_name() << std::endl;
			// exit(EXIT_FAILURE);

			// Keep at same point
			_vals[timepoint+1] = _vals[timepoint];

		} else {
			_vals[timepoint+1] = _vals[timepoint] + dt * _diff_eq->get_val_at_timepoint(timepoint);
		};
	};

	/********************
	Moment
	********************/

	std::shared_ptr<Moment> IxnParam::Impl::get_moment() const {
		return _moment;
	};

	/********************
	Adjoint
	********************/

	void IxnParam::Impl::set_adjoint(std::shared_ptr<Adjoint> adjoint) {
		_adjoint = adjoint;
	};
	std::shared_ptr<Adjoint> IxnParam::Impl::get_adjoint() const {
		return _adjoint;
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
		for (auto timepoint=0; timepoint<_no_timepoints; timepoint++) {
			f << _vals[timepoint] << "\n";
		};

		// Close
		f.close();
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
	Diff eq rhs that this ixn param appears in
	********************/

	void IxnParam::add_diff_eq_dependency(std::shared_ptr<DiffEqRHS> diff_eq, int idx) {
		_impl->add_diff_eq_dependency(diff_eq,idx);
	};
	const std::vector<std::pair<std::shared_ptr<DiffEqRHS>,int>>& IxnParam::get_diff_eq_dependencies() const {
		return _impl->get_diff_eq_dependencies();
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
	Fixed value to IC
	********************/

	void IxnParam::set_is_val_fixed_to_init_cond(bool fixed) {
		_impl->set_is_val_fixed_to_init_cond(fixed);
	};
	bool IxnParam::get_is_val_fixed_to_init_cond() const {
		return _impl->get_is_val_fixed_to_init_cond();
	};

	void IxnParam::set_are_vals_fixed(bool fixed) {
		_impl->set_are_vals_fixed(fixed);
	};
	bool IxnParam::get_are_vals_fixed() const {
		return _impl->get_are_vals_fixed();
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

	void IxnParam::set_val_at_timepoint(int timepoint, double val) {
		_impl->set_val_at_timepoint(timepoint,val);
	};
	double IxnParam::get_val_at_timepoint(int timepoint) const {
		return _impl->get_val_at_timepoint(timepoint);
	};

	/********************
	Diff eq
	********************/

	std::shared_ptr<DiffEqRHS> IxnParam::get_diff_eq_rhs() const {
		return _impl->get_diff_eq_rhs();
	};
	void IxnParam::set_diff_eq_rhs(std::shared_ptr<DiffEqRHS> diff_eq) {
		_impl->set_diff_eq_rhs(diff_eq);
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

	/********************
	Adjoint
	********************/

	void IxnParam::set_adjoint(std::shared_ptr<Adjoint> adjoint) {
		_impl->set_adjoint(adjoint);
	};
	std::shared_ptr<Adjoint> IxnParam::get_adjoint() const {
		return _impl->get_adjoint();
	};

	/********************
	Write to file
	********************/

	void IxnParam::write_to_file(std::string fname) const {
		_impl->write_to_file(fname);
	};

};





