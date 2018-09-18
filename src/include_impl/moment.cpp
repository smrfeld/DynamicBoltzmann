#include "../../include/dynamicboltz_bits/moment.hpp"

// Other headers
#include "../../include/dynamicboltz_bits/general.hpp"
#include "../../include/dynamicboltz_bits/ixn_param.hpp"
#include "../../include/dynamicboltz_bits/unit_visible.hpp"
#include "../../include/dynamicboltz_bits/connections.hpp"
#include "../../include/dynamicboltz_bits/species.hpp"

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

	Moment::Moment(IxnParamType type) {
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

		_type = other._type;

		_monitor_h = other._monitor_h; 
		_monitor_j = other._monitor_j; 
		_monitor_k = other._monitor_k; 

		_no_timesteps = other._no_timesteps;
		_no_timepoints = other._no_timepoints;
		_batch_size = other._batch_size;

		// reaped
		_vals_awake_reaped = other._vals_awake_reaped;
		_vals_asleep_reaped = other._vals_asleep_reaped;

		// averaged
		_vals_awake_averaged = other._vals_awake_averaged;
		_vals_asleep_averaged = other._vals_asleep_averaged;

		_sp_h = other._sp_h;
		_sp_b = other._sp_b;
		_sp_j = other._sp_j;
		_sp_k = other._sp_k;
		_sp_w = other._sp_w;

		// Reset the other

		other._monitor_h.clear(); 
		other._monitor_j.clear(); 
		other._monitor_k.clear(); 

		other._no_timesteps = 0;
		other._no_timepoints = 0;
		other._batch_size = 0;

		other._vals_awake_reaped = nullptr;
		other._vals_asleep_reaped = nullptr;

		other._vals_awake_averaged = nullptr;
		other._vals_asleep_averaged = nullptr;

		other._sp_h.clear();
		other._sp_b.clear();
		other._sp_j.clear();
		other._sp_k.clear();
		other._sp_w.clear();
	};
	void Moment::_copy(const Moment& other) {

		_type = other._type;

		_monitor_h = other._monitor_h; 
		_monitor_j = other._monitor_j; 
		_monitor_k = other._monitor_k; 

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

		_sp_h = other._sp_h;
		_sp_b = other._sp_b;
		_sp_j = other._sp_j;
		_sp_k = other._sp_k;
		_sp_w = other._sp_w;
	};

	/********************
	Check setup
	********************/

	void Moment::check_setup() const {
		// ...
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
	void Moment::add_conn_to_monitor_j(ConnVV *conn, bool this_order) {
		if (_type != IxnParamType::J) {
			std::cerr << ">>> Error: Moment::add_conn_to_monitor_j <<< tried to add connVV but moment is not of type J" << std::endl;
			exit(EXIT_FAILURE);
		};
		_monitor_j.push_back(std::make_pair(conn,this_order));
	};
	void Moment::add_conn_to_monitor_k(ConnVVV *conn, bool this_order) {
		if (_type != IxnParamType::K) {
			std::cerr << ">>> Error: Moment::add_conn_to_monitor_k <<< tried to add connVVV but moment is not of type K" << std::endl;
			exit(EXIT_FAILURE);
		};
		_monitor_k.push_back(std::make_pair(conn,this_order));
	};

	// Add species
	void Moment::add_species_h(Sptr species) {
		if (_type != IxnParamType::H) {
			std::cerr << ">>> ERROR: Moment::Impl::add_species_h <<< Not of type H." << std::endl;
			exit(EXIT_FAILURE);
		};
		_sp_h.push_back(species);	
	};
	void Moment::add_species_b(Sptr species) {
		if (_type != IxnParamType::B) {
			std::cerr << ">>> ERROR: Moment::Impl::add_species_h <<< Not of type B." << std::endl;
			exit(EXIT_FAILURE);
		};
		_sp_b.push_back(species);
	};
	void Moment::add_species_j(Sptr species1, Sptr species2) {
		if (_type != IxnParamType::J) {
			std::cerr << ">>> ERROR: Moment::Impl::add_species_h <<< Not of type J." << std::endl;
			exit(EXIT_FAILURE);
		};
		_sp_j.push_back(Sptr2(species1,species2));
	};
	void Moment::add_species_k(Sptr species1, Sptr species2, Sptr species3) {
		if (_type != IxnParamType::K) {
			std::cerr << ">>> ERROR: Moment::Impl::add_species_h <<< Not of type K." << std::endl;
			exit(EXIT_FAILURE);
		};
		_sp_k.push_back(Sptr3(species1,species2,species3));
	};
	void Moment::add_species_w(Sptr speciesV, Sptr speciesH) {
		if (_type != IxnParamType::W) {
			std::cerr << ">>> ERROR: Moment::Impl::add_species_h <<< Not of type W." << std::endl;
			exit(EXIT_FAILURE);
		};
		_sp_w.push_back(Sptr2(speciesV,speciesH));
	};

	/********************
	Get species
	********************/

	const std::vector<Sptr>& Moment::get_species_h() const {
		if (_type != IxnParamType::H) {
			std::cerr << ">>> ERROR: Moment::Impl::get_species_h <<< Not of type H." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _sp_h;
	};
	const std::vector<Sptr>& Moment::get_species_b() const {
		if (_type != IxnParamType::B) {
			std::cerr << ">>> ERROR: Moment::Impl::get_species_h <<< Not of type B." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _sp_b;
	};
	const std::vector<Sptr2>& Moment::get_species_j() const {
		if (_type != IxnParamType::J) {
			std::cerr << ">>> ERROR: Moment::Impl::get_species_h <<< Not of type J." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _sp_j;
	};
	const std::vector<Sptr3>& Moment::get_species_k() const {
		if (_type != IxnParamType::K) {
			std::cerr << ">>> ERROR: Moment::Impl::get_species_h <<< Not of type K." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _sp_k;
	};
	const std::vector<Sptr2>& Moment::get_species_w() const {
		if (_type != IxnParamType::W) {
			std::cerr << ">>> ERROR: Moment::Impl::get_species_h <<< Not of type W." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _sp_w;
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

		if (timepoint >= _no_timepoints) {
			std::cerr << ">>> Error Moment::reap_as_timepoint <<< Max timepoint for moment is: " << _no_timepoints << " but tried: " << timepoint << std::endl;
			exit(EXIT_FAILURE);
		};
		if (i_batch >= _batch_size) {
			std::cerr << ">>> Error Moment::reap_as_timepoint <<< Batch size for moment is: " << _batch_size << " but tried: " << i_batch << std::endl;
			exit(EXIT_FAILURE);
		};

		// Get all counts
		double count = 0.0;
		if (_type == IxnParamType::H) {
			// H
			if (binary) {
				Sptr sp_unit;
				for (auto const &s: _monitor_h) {
					sp_unit = s->get_b_mode_species();
					for (auto const &sp: _sp_h) {
						if (sp_unit == sp) {
							count += 1.0;
						};
					};
				};
			} else {
				for (auto const &s: _monitor_h) {
					for (auto const &sp: _sp_h) {
						count += s->get_p_mode_prob(sp);
					};
				};
			};
		} else if (_type == IxnParamType::J) {
			// J
			for (auto const &conn: _monitor_j) {
				for (auto &sp_pair: _sp_j) {
					count += conn.first->get_count(sp_pair.s1,sp_pair.s2,binary,conn.second);
				};
			};
		} else if (_type == IxnParamType::K) {
			// K
			for (auto const &conn: _monitor_k) {
				for (auto &sp_triplet: _sp_k) {
					count += conn.first->get_count(sp_triplet.s1,sp_triplet.s2,sp_triplet.s3,binary,conn.second);
				};
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


};


