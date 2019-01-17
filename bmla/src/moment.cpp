#include "../include/bmla_bits/moment.hpp"

// Other headers
#include "../include/bmla_bits/general.hpp"
#include "../include/bmla_bits/ixn_param.hpp"
#include "../include/bmla_bits/unit.hpp"
#include "../include/bmla_bits/connections.hpp"
#include "../include/bmla_bits/species.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

/************************************
* Namespace for bmla
************************************/

namespace bmla {

	/****************************************
	Moment
	****************************************/
	
	/********************
	Constructor
	********************/

	Moment::Moment(std::string name, IxnParamType type) {
		_name = name;
		_type = type;

		_batch_size = 1;
        _no_markov_chains = 1;
        
		// reaped
		_vals_awake_reaped = new double[_batch_size];
		_vals_asleep_reaped = new double[_no_markov_chains];
		std::fill_n(_vals_awake_reaped, _batch_size, 0.0);
		std::fill_n(_vals_asleep_reaped, _no_markov_chains, 0.0);

		// averaged
		_val_awake_averaged = 0.0;
		_val_asleep_averaged = 0.0;

		// Fixed awake moment
		_is_awake_moment_fixed = false;
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

		_batch_size = other._batch_size;
        _no_markov_chains = other._no_markov_chains;
        
		// reaped
		_vals_awake_reaped = other._vals_awake_reaped;
		_vals_asleep_reaped = other._vals_asleep_reaped;

		// averaged
		_val_awake_averaged = other._val_awake_averaged;
		_val_asleep_averaged = other._val_asleep_averaged;
		
		_is_awake_moment_fixed = other._is_awake_moment_fixed;

		// Reset the other
		other._name = "";

		other._monitor_h.clear(); 
		other._monitor_b.clear();
		other._monitor_j.clear(); 
		other._monitor_k.clear(); 
		other._monitor_w.clear(); 
		other._monitor_x.clear();

		other._batch_size = 0;
        other._no_markov_chains = 0;
        
		other._vals_awake_reaped = nullptr;
		other._vals_asleep_reaped = nullptr;

		other._val_awake_averaged = 0.0;
		other._val_asleep_averaged = 0.0;

		other._is_awake_moment_fixed = false;
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

		_batch_size = other._batch_size;
        _no_markov_chains = other._no_markov_chains;
        
		// reaped
		_vals_awake_reaped = new double[_batch_size];
		_vals_asleep_reaped = new double[_no_markov_chains];
		std::copy(other._vals_awake_reaped, other._vals_awake_reaped+_batch_size,_vals_awake_reaped);
		std::copy(other._vals_asleep_reaped, other._vals_asleep_reaped+_no_markov_chains,_vals_asleep_reaped);

		// averaged
		_val_awake_averaged = other._val_awake_averaged;
		_val_asleep_averaged = other._val_asleep_averaged;

		_is_awake_moment_fixed = other._is_awake_moment_fixed;
	};

	/********************
	Verbose
	********************/

	void Moment::print_moment_comparison() const {
		std::cout << "(" << _val_awake_averaged << "," << _val_asleep_averaged << ") " << std::endl;
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
	Batch size
	********************/

	int Moment::get_batch_size() const {
		return _batch_size;
	};
	void Moment::set_batch_size(int batch_size) {
		// Clear old
		
		// reaped
		safeDelArr(_vals_awake_reaped);

		// New
		
		// Params
		_batch_size = batch_size;

		// Reaped vals
		_vals_awake_reaped = new double[_batch_size];
		std::fill_n(_vals_awake_reaped,_batch_size,0.0);

		// Averaged vals
		_val_awake_averaged = 0.0;
	};

    int Moment::get_no_markov_chains() const {
        return _no_markov_chains;
    };
    void Moment::set_no_markov_chains(int no_markov_chains) {
        // Clear old
        
        // reaped
        safeDelArr(_vals_asleep_reaped);
        
        // New
        
        // Params
        _no_markov_chains = no_markov_chains;
        
        // Reaped vals
        _vals_asleep_reaped = new double[_batch_size];
        std::fill_n(_vals_asleep_reaped,_batch_size,0.0);
        
        // Averaged vals
        _val_asleep_averaged = 0.0;
    };

    
	/********************
	Reset
	********************/

	void Moment::reset_to_zero(MomentType type) {
		if (type == MomentType::AWAKE) {
			for (auto i_batch=0; i_batch<_batch_size; i_batch++) {
				set_moment_sample(MomentType::AWAKE, i_batch, 0.0);
			};
			set_moment(MomentType::AWAKE, 0.0);
		} else {
			for (auto i_chain=0; i_chain<_no_markov_chains; i_chain++) {
				set_moment_sample(MomentType::ASLEEP, i_chain, 0.0);
			};
			set_moment(MomentType::ASLEEP, 0.0);
		};
	};

	/********************
	Fixed awake
	********************/

	void Moment::set_is_awake_moment_fixed(bool flag) {
		_is_awake_moment_fixed = flag;
	};
	bool Moment::get_is_awake_moment_fixed() const {
		return _is_awake_moment_fixed;
	};

	/********************
	Get/set moment
	********************/

	double Moment::get_moment(MomentType type) const {
		if (type == MomentType::AWAKE) {
			return _val_awake_averaged;
		} else {
			return _val_asleep_averaged;
		};
	};
	void Moment::set_moment(MomentType type, double val) {
		if (type == MomentType::AWAKE) {
			_val_awake_averaged = val;
		} else {
			_val_asleep_averaged = val;
		};
	};

	// Batch
	double Moment::get_moment_sample(MomentType type, int i_sample) const {
		if (type == MomentType::AWAKE) {
			return _vals_awake_reaped[i_sample];
		} else {
			return _vals_asleep_reaped[i_sample];
		};	
	};
	void Moment::set_moment_sample(MomentType type, int i_sample, double val) {
		if (type == MomentType::AWAKE) {
			_vals_awake_reaped[i_sample] = val;
		} else {
			_vals_asleep_reaped[i_sample] = val;
		};
	};

	/********************
	Reap from sampler
	********************/

	void Moment::reap_sample(MomentType type, int i_sample) {

        if (type == MomentType::AWAKE) {
            if (i_sample >= _batch_size) {
                std::cerr << ">>> Error Moment::reap_in_batch <<< Batch size for awake moment is: " << _batch_size << " but tried: " << i_sample << std::endl;
                exit(EXIT_FAILURE);
            };
        } else {
            if (i_sample >= _no_markov_chains) {
                std::cerr << ">>> Error Moment::reap_in_batch <<< No markov chains for  asleep moment is: " << _no_markov_chains << " but tried: " << i_sample << std::endl;
                exit(EXIT_FAILURE);
            };
        };

		// Get all counts
		double count = 0.0;
		if (_type == IxnParamType::H) {
			// H
			for (auto const &unit: _monitor_h) {
				count += unit->get_moment(_name);
			};
		} else if (_type == IxnParamType::B) {
			// B
			for (auto const &unit: _monitor_b) {
				count += unit->get_moment(_name);
			};
		} else if (_type == IxnParamType::J) {
			// J
			for (auto const &conn: _monitor_j) {
				count += conn->get_moment(_name);
			};
		} else if (_type == IxnParamType::K) {
			// K
			for (auto const &conn: _monitor_k) {
				count += conn->get_moment(_name);
			};
		} else if (_type == IxnParamType::W) {
			// K
			for (auto const &conn: _monitor_w) {
				count += conn->get_moment(_name);
			};
		} else if (_type == IxnParamType::X) {
			// K
			for (auto const &conn: _monitor_x) {
				count += conn->get_moment(_name);
			};
		};

		// Set
		set_moment_sample(type,i_sample,count);
	};

	/********************
	Average reaps
	********************/

	void Moment::average_reaps(MomentType type) {
		if (type == MomentType::AWAKE) {

			// Awake moment
			_val_awake_averaged = 0.;
			for (auto i=0; i<_batch_size; i++) {
				_val_awake_averaged += _vals_awake_reaped[i];
			};
			_val_awake_averaged /= _batch_size;

		} else {

			// Asleep moment
			_val_asleep_averaged = 0.;
			for (auto i=0; i<_no_markov_chains; i++) {
				_val_asleep_averaged += _vals_asleep_reaped[i];
			};
			_val_asleep_averaged /= _no_markov_chains;

		};
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
		f << _val_awake_averaged << " " << _val_asleep_averaged << "\n";

		// Close
		f.close();
	};

};


