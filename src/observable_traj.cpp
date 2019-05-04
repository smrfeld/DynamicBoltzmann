#include "../include/dblz_bits/observable_traj.hpp"

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
	ObservableTraj
	****************************************/
	
	/********************
	Constructor
	********************/

	ObservableTraj::ObservableTraj(ITptr ixn_param, int layer, Sptr species) {
        _ixn_param = ixn_param;
        _layer = layer;
        _species = species;
        
        set_no_timesteps(0);
    };
	ObservableTraj::ObservableTraj(const ObservableTraj& other) {
		_copy(other);
	};
	ObservableTraj::ObservableTraj(ObservableTraj&& other) {
		_move(other);
	};
	ObservableTraj& ObservableTraj::operator=(const ObservableTraj& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	ObservableTraj& ObservableTraj::operator=(ObservableTraj&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
	};
	ObservableTraj::~ObservableTraj() {
		_clean_up();
	};

	void ObservableTraj::_clean_up() {
    };
	void ObservableTraj::_move(ObservableTraj &other) {
        _ixn_param = other._ixn_param;
        _layer = other._layer;
        _species = other._species;
        _no_timesteps = other._no_timesteps;
        _no_timepoints = other._no_timepoints;
        _vals = other._vals;
        
        other._ixn_param = nullptr;
        other._layer = 0;
        other._species = nullptr;
        other._no_timesteps = 0;
        other._no_timepoints = 0;
        other._vals.clear();
	};
	void ObservableTraj::_copy(const ObservableTraj& other) {
        _ixn_param = other._ixn_param;
        _layer = other._layer;
        _species = other._species;
        _no_timesteps = other._no_timesteps;
        _no_timepoints = other._no_timepoints;
        _vals = other._vals;
	};

	/********************
	Name
	********************/

    ITptr ObservableTraj::get_ixn_param_traj() const {
        return _ixn_param;
    };
    int ObservableTraj::get_layer() const {
        return _layer;
    };
    Sptr ObservableTraj::get_species() const {
        return _species;
    };

    /********************
     Timepoints
     ********************/
    
    void ObservableTraj::set_no_timesteps(int no_timesteps) {
        _no_timesteps = no_timesteps;
        _no_timepoints = no_timesteps + 1;
        
        while (_vals.size() < _no_timepoints) {
            _vals.push_back(0.0);
        };
        while (_vals.size() > _no_timepoints) {
            _vals.pop_back();
        };
    };
    int ObservableTraj::get_no_timesteps() const {
        return _no_timesteps;
    };

	/********************
	Get/set moment
	********************/

    // Get moment
    double ObservableTraj::get_val_at_timepoint(int timepoint) const {
        return _vals.at(timepoint);
    };
    void ObservableTraj::increment_val_at_timepoint(int timepoint, double val) {
        _vals[timepoint] += val;
    };
    void ObservableTraj::set_val_at_timepoint(int timepoint, double val) {
        _vals[timepoint] = val;
    };
    void ObservableTraj::reset_val_at_timepoint(int timepoint) {
        _vals[timepoint] = 0.0;
    };
};


