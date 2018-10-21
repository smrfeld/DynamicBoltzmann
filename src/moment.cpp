#include "../include/dynamicboltz_bits/moment.hpp"

// Other headers
#include "../include/dynamicboltz_bits/general.hpp"
#include "../include/dynamicboltz_bits/ixn_param.hpp"
#include "../include/dynamicboltz_bits/unit.hpp"
#include "../include/dynamicboltz_bits/connections.hpp"
#include "../include/dynamicboltz_bits/species.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

/************************************
* Namespace for dblz
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

		_no_timesteps = 0;
		_no_timepoints = 1;

		_batch_size = 1;

		// reaped
		_vals_awake_reaped = new double[_no_timepoints*_batch_size];
		_vals_asleep_reaped = new double[_no_timepoints*_batch_size];
		std::fill_n(_vals_awake_reaped, _batch_size*_no_timepoints, 0.0);
		std::fill_n(_vals_asleep_reaped, _batch_size*_no_timepoints, 0.0);

		// averaged
		_vals_awake_averaged = new double[_no_timepoints];
		_vals_asleep_averaged = new double[_no_timepoints];
		std::fill_n(_vals_awake_averaged, _no_timepoints, 0.0);
		std::fill_n(_vals_asleep_averaged, _no_timepoints, 0.0);
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

		_monitor_h = other._monitor_h; 
		_monitor_b = other._monitor_b; 
		_monitor_j = other._monitor_j; 
		_monitor_k = other._monitor_k; 
		_monitor_w = other._monitor_w; 
		_monitor_x = other._monitor_x; 

		_no_timesteps = other._no_timesteps;
		_no_timepoints = other._no_timepoints;
		_batch_size = other._batch_size;

		// reaped
		_vals_awake_reaped = other._vals_awake_reaped;
		_vals_asleep_reaped = other._vals_asleep_reaped;

		// averaged
		_vals_awake_averaged = other._vals_awake_averaged;
		_vals_asleep_averaged = other._vals_asleep_averaged;
		
		// Reset the other
		other._name = "";

		other._monitor_h.clear(); 
		other._monitor_b.clear();
		other._monitor_j.clear(); 
		other._monitor_k.clear(); 
		other._monitor_w.clear(); 
		other._monitor_x.clear(); 

		other._no_timesteps = 0;
		other._no_timepoints = 0;
		other._batch_size = 0;

		other._vals_awake_reaped = nullptr;
		other._vals_asleep_reaped = nullptr;

		other._vals_awake_averaged = nullptr;
		other._vals_asleep_averaged = nullptr;
	};
	void Moment::_copy(const Moment& other) {
		_name = other._name;
		_type = other._type;

		_monitor_h = other._monitor_h; 
		_monitor_b = other._monitor_b;
		_monitor_j = other._monitor_j; 
		_monitor_k = other._monitor_k; 
		_monitor_w = other._monitor_w; 
		_monitor_x = other._monitor_x; 

		_no_timesteps = other._no_timesteps;
		_no_timepoints = other._no_timepoints;
		_batch_size = other._batch_size;

		// reaped
		_vals_awake_reaped = new double[_no_timepoints*_batch_size];
		_vals_asleep_reaped = new double[_no_timepoints*_batch_size];
		std::copy(other._vals_awake_reaped, other._vals_awake_reaped+_no_timepoints*_batch_size,_vals_awake_reaped);
		std::copy(other._vals_asleep_reaped, other._vals_asleep_reaped+_no_timepoints*_batch_size,_vals_asleep_reaped);

		// averaged
		_vals_awake_averaged = new double[_no_timepoints];
		std::copy( other._vals_awake_averaged, other._vals_awake_averaged + _no_timepoints, _vals_awake_averaged );
		_vals_asleep_averaged = new double[_no_timepoints];
		std::copy( other._vals_asleep_averaged, other._vals_asleep_averaged + _no_timepoints, _vals_asleep_averaged );
	};

	/********************
	Verbose
	********************/

	void Moment::print_moment_comparison() const {
		// eg timesteps = 2; timepoints = 3; 0->1, 1->2; three loops
		for (auto timepoint=0; timepoint<=_no_timepoints-1; timepoint++) {
			std::cout << "(" << _vals_awake_averaged[timepoint] << "," << _vals_asleep_averaged[timepoint] << ") ";
		};
		std::cout << std::endl;
	};

	/********************
	Finish setup
	********************/

	void Moment::add_unit_to_monitor_h(UnitVisible *uv) {
		if (_type != IxnParamType::H) {
			std::cerr << ">>> Error: Moment::add_unit_to_monitor_h <<< tried to add visible unit but moment is not of type H" << std::endl;
			exit(EXIT_FAILURE);
		};
		_monitor_h.push_back(uv);
	};
	void Moment::add_unit_to_monitor_b(UnitHidden *uh) {
		if (_type != IxnParamType::B) {
			std::cerr << ">>> Error: Moment::add_unit_to_monitor_b <<< tried to add visible unit but moment is not of type B" << std::endl;
			exit(EXIT_FAILURE);
		};
		_monitor_b.push_back(uh);
	};
	void Moment::add_conn_to_monitor_j(ConnVV *conn) {
		if (_type != IxnParamType::J) {
			std::cerr << ">>> Error: Moment::add_conn_to_monitor_j <<< tried to add connVV but moment is not of type J" << std::endl;
			exit(EXIT_FAILURE);
		};
		_monitor_j.push_back(conn);
	};
	void Moment::add_conn_to_monitor_k(ConnVVV *conn) {
		if (_type != IxnParamType::K) {
			std::cerr << ">>> Error: Moment::add_conn_to_monitor_k <<< tried to add connVVV but moment is not of type K" << std::endl;
			exit(EXIT_FAILURE);
		};
		_monitor_k.push_back(conn);
	};
	void Moment::add_conn_to_monitor_w(ConnVH *conn) {
		if (_type != IxnParamType::W) {
			std::cerr << ">>> Error: Moment::add_conn_to_monitor_w <<< tried to add connVH but moment is not of type W" << std::endl;
			exit(EXIT_FAILURE);
		};
		_monitor_w.push_back(conn);
	};
	void Moment::add_conn_to_monitor_x(ConnHH *conn) {
		if (_type != IxnParamType::X) {
			std::cerr << ">>> Error: Moment::add_conn_to_monitor_x <<< tried to add connHH but moment is not of type X" << std::endl;
			exit(EXIT_FAILURE);
		};
		_monitor_x.push_back(conn);
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
	Number timesteps
	********************/

	int Moment::get_no_timesteps() const {
		return _no_timesteps;
	};
	void Moment::set_no_timesteps(int no_timesteps) {
		// Clear old
		
		// reaped
		safeDelArr(_vals_awake_reaped);
		safeDelArr(_vals_asleep_reaped);

		// averaged
		safeDelArr(_vals_awake_averaged);
		safeDelArr(_vals_asleep_averaged);

		// New
		
		// Params
		_no_timesteps = no_timesteps;
		_no_timepoints = _no_timesteps+1;

		// Reaped vals
		_vals_awake_reaped = new double[_no_timepoints*_batch_size];
		_vals_asleep_reaped = new double[_no_timepoints*_batch_size];
		std::fill_n(_vals_awake_reaped,_no_timepoints*_batch_size,0.0);
		std::fill_n(_vals_asleep_reaped,_no_timepoints*_batch_size,0.0);

		// Averaged vals
		_vals_awake_averaged = new double[_no_timepoints];
		_vals_asleep_averaged = new double[_no_timepoints];
		std::fill_n(_vals_awake_averaged,_no_timepoints,0.0);
		std::fill_n(_vals_asleep_averaged,_no_timepoints,0.0);
	};

	/********************
	Batch size
	********************/

	int Moment::get_batch_size() const {
		return _batch_size;
	};
	void Moment::set_batch_size(int batch_size) {
		// Clear old
		
		// reaped
		safeDelArr(_vals_awake_reaped);
		safeDelArr(_vals_asleep_reaped);

		// averaged
		safeDelArr(_vals_awake_averaged);
		safeDelArr(_vals_asleep_averaged);

		// New
		
		// Params
		_batch_size = batch_size;

		// Reaped vals
		_vals_awake_reaped = new double[_no_timepoints*_batch_size];
		_vals_asleep_reaped = new double[_no_timepoints*_batch_size];
		std::fill_n(_vals_awake_reaped,_no_timepoints*_batch_size,0.0);
		std::fill_n(_vals_asleep_reaped,_no_timepoints*_batch_size,0.0);

		// Averaged vals
		_vals_awake_averaged = new double[_no_timepoints];
		_vals_asleep_averaged = new double[_no_timepoints];
		std::fill_n(_vals_awake_averaged,_no_timepoints,0.0);
		std::fill_n(_vals_asleep_averaged,_no_timepoints,0.0);
	};

	/********************
	Reset
	********************/

	void Moment::reset_to_zero() {
		for (auto timepoint=0; timepoint<_no_timepoints; timepoint++) {
			for (auto i_batch=0; i_batch<_batch_size; i_batch++) {
				set_moment_at_timepoint_in_batch(MomentType::AWAKE, timepoint, i_batch, 0.0);
				set_moment_at_timepoint_in_batch(MomentType::ASLEEP, timepoint, i_batch, 0.0);
			};
			set_moment_at_timepoint(MomentType::AWAKE, timepoint, 0.0);
			set_moment_at_timepoint(MomentType::ASLEEP, timepoint, 0.0);
		};
	};

	/********************
	Get/set moment
	********************/

	double Moment::get_moment_at_timepoint(MomentType type, int timepoint) const {
		if (type == MomentType::AWAKE) {
			return _vals_awake_averaged[timepoint];
		} else {
			return _vals_asleep_averaged[timepoint];
		};
	};
	void Moment::set_moment_at_timepoint(MomentType type, int timepoint, double val) {
		if (type == MomentType::AWAKE) {
			_vals_awake_averaged[timepoint] = val;
		} else {
			_vals_asleep_averaged[timepoint] = val;
		};
	};

	// Batch
	double Moment::get_moment_at_timepoint_in_batch(MomentType type, int timepoint, int i_batch) const {
		if (type == MomentType::AWAKE) {
			return _vals_awake_reaped[timepoint*_batch_size + i_batch];
		} else {
			return _vals_asleep_reaped[timepoint*_batch_size + i_batch];
		};	
	};
	void Moment::set_moment_at_timepoint_in_batch(MomentType type, int timepoint, int i_batch, double val) {
		if (type == MomentType::AWAKE) {
			_vals_awake_reaped[timepoint*_batch_size + i_batch] = val;
		} else {
			_vals_asleep_reaped[timepoint*_batch_size + i_batch] = val;
		};
	};

	/********************
	Reap from sampler
	********************/

	void Moment::reap_as_timepoint_in_batch(MomentType type, int timepoint, int i_batch, bool binary) {

		/*
		if (timepoint >= _no_timepoints) {
			std::cerr << ">>> Error Moment::reap_as_timepoint <<< Max timepoint for moment is: " << _no_timepoints << " but tried: " << timepoint << std::endl;
			exit(EXIT_FAILURE);
		};
		if (i_batch >= _batch_size) {
			std::cerr << ">>> Error Moment::reap_as_timepoint <<< Batch size for moment is: " << _batch_size << " but tried: " << i_batch << std::endl;
			exit(EXIT_FAILURE);
		};
		*/

		// Get all counts
		double count = 0.0;
		if (_type == IxnParamType::H) {
			// H
			for (auto const &unit: _monitor_h) {
				count += unit->get_moment(_name,binary);
			};
		} else if (_type == IxnParamType::B) {
			// B
			for (auto const &unit: _monitor_b) {
				count += unit->get_moment(_name,binary);
			};
		} else if (_type == IxnParamType::J) {
			// J
			for (auto const &conn: _monitor_j) {
				count += conn->get_moment(_name,binary);
			};
		} else if (_type == IxnParamType::K) {
			// K
			for (auto const &conn: _monitor_k) {
				count += conn->get_moment(_name,binary);
			};
		} else if (_type == IxnParamType::W) {
			// W
			for (auto const &conn: _monitor_w) {
				count += conn->get_moment(_name,binary);
			};
		} else if (_type == IxnParamType::X) {
			// X
			for (auto const &conn: _monitor_x) {
				count += conn->get_moment(_name,binary);
			};
		};

		// Set
		set_moment_at_timepoint_in_batch(type,timepoint,i_batch,count);
	};

	/********************
	Average reaps
	********************/

	void Moment::average_reaps_as_timepoint(MomentType type, int timepoint) {
		if (type == MomentType::AWAKE) {

			// Awake moment
			_vals_awake_averaged[timepoint] = 0.;
			for (auto i=0; i<_batch_size; i++) {
				_vals_awake_averaged[timepoint] += _vals_awake_reaped[timepoint*_batch_size + i];
			};
			_vals_awake_averaged[timepoint] /= _batch_size;

		} else {

			// Asleep moment
			_vals_asleep_averaged[timepoint] = 0.;
			for (auto i=0; i<_batch_size; i++) {
				_vals_asleep_averaged[timepoint] += _vals_asleep_reaped[timepoint*_batch_size + i];
			};
			_vals_asleep_averaged[timepoint] /= _batch_size;

		};
	};

	/********************
	Write
	********************/

	void Moment::write_to_file(std::string fname) const {
		std::ofstream f;

		// Open
		f.open(fname);

		// Make sure we found it
		if (!f.is_open()) {
			std::cerr << ">>> Error: Moment::write_to_file <<< could not write to file: " << fname << std::endl;
			exit(EXIT_FAILURE);
		};

		// Go through all time
		for (auto timepoint=0; timepoint<_no_timepoints; timepoint++) {
			f << _vals_awake_averaged[timepoint] << " " << _vals_asleep_averaged[timepoint] << "\n";
		};

		// Close
		f.close();
	};

};


