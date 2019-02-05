#include "../include/dblz_bits/lattice_traj.hpp"

// Other headers
#include "../include/dblz_bits/general.hpp"
#include "../include/dblz_bits/species.hpp"
#include "../include/dblz_bits/moment.hpp"
#include "../include/dblz_bits/ixn_param.hpp"
#include "../include/dblz_bits/ixn_param_traj.hpp"
#include "../include/dblz_bits/lattice.hpp"

#include <iostream>
#include <fstream>
#include <numeric>
#include "math.h"
#include <ctime>
#include <sstream>
#include <random>
#include <algorithm>

/************************************
* Namespace for bmla
************************************/

namespace dblz {

	/****************************************
	LatticeTraj
	****************************************/

	/********************
	Constructor
	********************/

	LatticeTraj::LatticeTraj(int no_dims, int box_length, std::vector<Sptr> species_visible)
	{
        // Make first lattice
        _lattices.push_back(std::make_shared<Lattice>(no_dims,box_length,species_visible,false));
        
        // Set no timesteps/timepoints
        set_no_timesteps(0);
	};
	LatticeTraj::LatticeTraj(const LatticeTraj& other) {
		_copy(other);
	};
	LatticeTraj::LatticeTraj(LatticeTraj&& other) {
		_move(other);
	};
	LatticeTraj& LatticeTraj::operator=(const LatticeTraj& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	LatticeTraj& LatticeTraj::operator=(LatticeTraj&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
	};
	LatticeTraj::~LatticeTraj() {
		_clean_up();
	};

	void LatticeTraj::_clean_up() {
        // Nothing...
	};
	void LatticeTraj::_copy(const LatticeTraj& other) {
        for (auto l: other._lattices) {
            _lattices.push_back(std::make_shared<Lattice>(*l));
        };
        _no_timepoints = other._no_timepoints;
        _no_timesteps = other._no_timesteps;
        _ixn_param_trajs = other._ixn_param_trajs;
    };
	void LatticeTraj::_move(LatticeTraj& other) {
        _lattices = other._lattices;
        _no_timepoints = other._no_timepoints;
        _no_timesteps = other._no_timesteps;
        _ixn_param_trajs = other._ixn_param_trajs;

		// Reset other
        other._lattices.clear();
        other._no_timesteps = 0;
        other._no_timepoints = 0;
        other._ixn_param_trajs.clear();
    };

    // ***************
    // MARK: - [Private] Add ixn param traj
    // ***************
    
    void LatticeTraj::_add_ixn_param_traj(ITptr ixn_param_traj) {
        auto it = std::find(_ixn_param_trajs.begin(),_ixn_param_trajs.end(),ixn_param_traj);
        if (it == _ixn_param_trajs.end()) {
            _ixn_param_trajs.push_back(ixn_param_traj);
        };
    };

    // ***************
    // MARK: - No timesteps
    // ***************
    
    int LatticeTraj::get_no_timesteps() const {
        return _no_timesteps;
    };
    void LatticeTraj::set_no_timesteps(int no_timesteps) {
        _no_timesteps = no_timesteps;
        _no_timepoints = _no_timesteps + 1;
        
        if (_no_timepoints <= 0) {
            std::cerr << ">>> LatticeTraj::set_no_timesteps <<< Error: no timepoints must be > 0! Now it is: " << _no_timepoints << std::endl;
            exit(EXIT_FAILURE);
        };
        
        while (_lattices.size() < _no_timepoints) {
            // Add
            _lattices.push_back(std::make_shared<Lattice>(*_lattices.back()));
        };
        while (_lattices.size() > _no_timepoints) {
            // Remove
            _lattices.pop_back();
        };
    };

    // ***************
    // MARK: - Getters
    // ***************
    
    int LatticeTraj::get_no_dims() const {
        return _lattices.front()->get_no_dims();
    };
    int LatticeTraj::get_box_length() const {
        return _lattices.front()->get_box_length();
    };
    
    int LatticeTraj::get_no_units_in_layer(int layer) const {
        return _lattices.front()->get_no_units_in_layer(layer);
    };
    
    int LatticeTraj::get_no_layers() const {
        return _lattices.front()->get_no_layers();
    };
    
    // Get all ixns
    const std::vector<ITptr>& LatticeTraj::get_all_ixn_param_trajs() const {
        return _ixn_param_trajs;
    };
    
    // ***************
    // MARK: - Markov chains
    // ***************
    
    int LatticeTraj::get_no_markov_chains(MCType type) const {
        return _lattices.front()->get_no_markov_chains(type);
    };

    void LatticeTraj::set_no_markov_chains(MCType type, int no_markov_chains) {
        for (auto l: _lattices) {
            l->set_no_markov_chains(type, no_markov_chains);
        };
    };
    
    // ***************
    // MARK: - Setup all timepoints
    // ***************
    
     // Add a layer
    void LatticeTraj::add_layer(int layer, int box_length, std::vector<Sptr> species) {
        for (auto l: _lattices) {
            l->add_layer(layer, box_length, species);
        };
    };
    
    // Biases
	void LatticeTraj::add_bias_all_layers(Sptr sp, ITptr bias) {
        _add_ixn_param_traj(bias);
        
        for (auto timepoint=0; timepoint<_no_timepoints; timepoint++) {
            _lattices.at(timepoint)->add_bias_all_layers(sp, bias->get_ixn_param_at_timepoint(timepoint));
        };
    };

    void LatticeTraj::add_bias_to_layer(int layer, Sptr sp, ITptr bias) {
        _add_ixn_param_traj(bias);

        for (auto timepoint=0; timepoint<_no_timepoints; timepoint++) {
            _lattices.at(timepoint)->add_bias_to_layer(layer, sp, bias->get_ixn_param_at_timepoint(timepoint));
        };
    };

    // Ixns
    void LatticeTraj::add_ixn_between_layers(int layer1, Sptr sp1, int layer2, Sptr sp2, ITptr ixn) {
        _add_ixn_param_traj(ixn);

        for (auto timepoint=0; timepoint<_no_timepoints; timepoint++) {
            _lattices.at(timepoint)->add_ixn_between_layers(layer1, sp1, layer2, sp2, ixn->get_ixn_param_at_timepoint(timepoint));
        };
    };
    
    // Set multiplier
    void LatticeTraj::set_multiplier_between_layers(int from_layer, int to_layer, double multiplier) {
        for (auto l: _lattices) {
            l->set_multiplier_between_layers(from_layer, to_layer, multiplier);
        };
    };
    void LatticeTraj::set_multiplier_for_bias_in_layer(int layer, double multiplier) {
        for (auto l: _lattices) {
            l->set_multiplier_for_bias_in_layer(layer, multiplier);
        };
    };
    
    // Add connections
    void LatticeTraj::add_conn(int layer1, int x1, int layer2, int x2) {
        for (auto l: _lattices) {
            l->add_conn(layer1,x1,layer2,x2);
        };
    };
    void LatticeTraj::add_conn(int layer1, int x1, int y1, int layer2, int x2, int y2) {
        for (auto l: _lattices) {
            l->add_conn(layer1,x1,y1,layer2,x2,y2);
        };
    };
    void LatticeTraj::add_conn(int layer1, int x1, int y1, int z1, int layer2, int x2, int y2, int z2) {
        for (auto l: _lattices) {
            l->add_conn(layer1,x1,y1,z1,layer2,x2,y2,z2);
        };
    };
};
