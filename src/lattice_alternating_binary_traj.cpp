#include "../include/dblz_bits/lattice_alternating_binary_traj.hpp"

// Other headers
#include "../include/dblz_bits/general.hpp"
#include "../include/dblz_bits/species.hpp"
#include "../include/dblz_bits/moment_diff.hpp"
#include "../include/dblz_bits/ixn_param.hpp"
#include "../include/dblz_bits/ixn_param_traj.hpp"
#include "../include/dblz_bits/lattice_alternating_binary.hpp"
#include "../include/dblz_bits/adjoint.hpp"

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
	LatticeAlternatingBinaryTraj
	****************************************/

	/********************
	Constructor
	********************/

    LatticeAlternatingBinaryTraj::LatticeAlternatingBinaryTraj(int no_dims, int box_length, std::vector<Sptr> species_visible)
	{
        _lattices[0] = std::make_shared<LatticeAlternatingBinary>(no_dims,box_length,species_visible);
        
        // Set no timesteps/timepoints
        set_no_timesteps(0,0);
	};
	LatticeAlternatingBinaryTraj::LatticeAlternatingBinaryTraj(const LatticeAlternatingBinaryTraj& other) {
		_copy(other);
	};
	LatticeAlternatingBinaryTraj::LatticeAlternatingBinaryTraj(LatticeAlternatingBinaryTraj&& other) {
		_move(other);
	};
	LatticeAlternatingBinaryTraj& LatticeAlternatingBinaryTraj::operator=(const LatticeAlternatingBinaryTraj& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	LatticeAlternatingBinaryTraj& LatticeAlternatingBinaryTraj::operator=(LatticeAlternatingBinaryTraj&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
	};
	LatticeAlternatingBinaryTraj::~LatticeAlternatingBinaryTraj() {
		_clean_up();
	};

	void LatticeAlternatingBinaryTraj::_clean_up() {
        // Nothing...
	};
	void LatticeAlternatingBinaryTraj::_copy(const LatticeAlternatingBinaryTraj& other) {
        for (auto l: other._lattices) {
            _lattices[l.first] = std::make_shared<LatticeAlternatingBinary>(*l.second);
        };
        _ixn_param_trajs = other._ixn_param_trajs;
        _bias_dict = other._bias_dict;
        _o2_ixn_dict = other._o2_ixn_dict;
    };
	void LatticeAlternatingBinaryTraj::_move(LatticeAlternatingBinaryTraj& other) {
        _lattices = other._lattices;
        _ixn_param_trajs = other._ixn_param_trajs;
        _bias_dict = other._bias_dict;
        _o2_ixn_dict = other._o2_ixn_dict;

		// Reset other
        other._lattices.clear();
        other._ixn_param_trajs.clear();
        other._bias_dict.clear();
        other._o2_ixn_dict.clear();
    };

    // ***************
    // MARK: - [Private] Add ixn param traj
    // ***************
    
    void LatticeAlternatingBinaryTraj::_add_ixn_param_traj(ITptr ixn_param_traj) {
        auto it = std::find(_ixn_param_trajs.begin(),_ixn_param_trajs.end(),ixn_param_traj);
        if (it == _ixn_param_trajs.end()) {
            _ixn_param_trajs.push_back(ixn_param_traj);
        };
    };

    // ***************
    // MARK: - No timesteps
    // ***************
    
    void LatticeAlternatingBinaryTraj::set_no_timesteps(int timepoint_start, int no_timesteps) {
        if (no_timesteps < 0) {
            std::cerr << ">>> LatticeAlternatingBinaryTraj::set_no_timesteps <<< Error: no timepoints must be > 0!" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        // Save at least one lattice
        auto itl = _lattices.end();
        itl--;
        auto lptr = itl->second;
        
        // Delete
        auto it = _lattices.begin();
        while (it != _lattices.end()) {
            if (it->first < timepoint_start || it->first > timepoint_start + no_timesteps) {
                // Delete
                it = _lattices.erase(it);
            } else {
                it++;
            };
        };
        
        // Add as needed
        for (auto timepoint=timepoint_start; timepoint<=timepoint_start + no_timesteps; timepoint++) {
            auto itf = _lattices.find(timepoint);
            if (itf == _lattices.end()) {
                // Does not exist; need to create
                _lattices[timepoint] = std::make_shared<LatticeAlternatingBinary>(*lptr);
                
                // Set correct ixns
                for (auto pr1: _bias_dict) {
                    for (auto pr2: pr1.second) {
                        _lattices[timepoint]->set_bias_of_layer(pr1.first, pr2.first, pr2.second->get_ixn_param_at_timepoint(timepoint));
                    };
                };
                for (auto pr1: _o2_ixn_dict) {
                    for (auto pr2: pr1.second) {
                        for (auto pr3: pr2.second) {
                            for (auto pr4: pr3.second) {
                                _lattices[timepoint]->set_ixn_between_layers(pr1.first, pr2.first, pr3.first, pr4.first, pr4.second->get_ixn_param_at_timepoint(timepoint));
                            };
                        };
                    };
                };
            };
        };
    };

    // ***************
    // MARK: - Getters
    // ***************
    
    int LatticeAlternatingBinaryTraj::get_no_dims() const {
        return _lattices.begin()->second->get_no_dims();
    };
    int LatticeAlternatingBinaryTraj::get_box_length() const {
        return _lattices.begin()->second->get_box_length();
    };
    
    int LatticeAlternatingBinaryTraj::get_no_units_in_layer(int layer) const {
        return _lattices.begin()->second->get_no_units_in_layer(layer);
    };
    
    int LatticeAlternatingBinaryTraj::get_no_layers() const {
        return _lattices.begin()->second->get_no_layers();
    };
    
    // Get all ixns
    const std::vector<ITptr>& LatticeAlternatingBinaryTraj::get_all_ixn_param_trajs() const {
        return _ixn_param_trajs;
    };
    
    // ***************
    // MARK: - Markov chains
    // ***************
    
    int LatticeAlternatingBinaryTraj::get_no_markov_chains(MCType type) const {
        return _lattices.begin()->second->get_no_markov_chains(type);
    };

    void LatticeAlternatingBinaryTraj::set_no_markov_chains(MCType type, int no_markov_chains) {
        for (auto l: _lattices) {
            l.second->set_no_markov_chains(type, no_markov_chains);
        };
    };
    
    // ***************
    // MARK: - Get lattice
    // ***************
    
    std::shared_ptr<LatticeAlternatingBinary> LatticeAlternatingBinaryTraj::get_lattice_at_timepoint(int timepoint) const {
        return _lattices.at(timepoint);
    };
    
    // ***************
    // MARK: - Setup all timepoints
    // ***************
    
     // Add a layer
    void LatticeAlternatingBinaryTraj::add_layer(int layer, int box_length, std::vector<Sptr> species) {
        for (auto l: _lattices) {
            l.second->add_layer(layer, box_length, species);
        };
    };

    // Biases
    void LatticeAlternatingBinaryTraj::set_bias_of_layer(int layer, Sptr sp, ITptr bias) {
        _bias_dict[layer][sp] = bias;

        _add_ixn_param_traj(bias);

        for (auto l: _lattices) {
            l.second->set_bias_of_layer(layer, sp, bias->get_ixn_param_at_timepoint(l.first));
        };
    };

    // Ixns
    void LatticeAlternatingBinaryTraj::set_ixn_between_layers(int layer1, Sptr sp1, int layer2, Sptr sp2, ITptr ixn) {
        _o2_ixn_dict[layer1][sp1][layer2][sp2] = ixn;
        _o2_ixn_dict[layer2][sp2][layer1][sp1] = ixn;

        _add_ixn_param_traj(ixn);

        for (auto l: _lattices) {
            l.second->set_ixn_between_layers(layer1, sp1, layer2, sp2, ixn->get_ixn_param_at_timepoint(l.first));
        };
    };
    
    // Set multiplier
    void LatticeAlternatingBinaryTraj::set_multiplier_between_layers(int from_layer, int to_layer, double multiplier) {
        for (auto l: _lattices) {
            l.second->set_multiplier_between_layers(from_layer, to_layer, multiplier);
        };
    };
    void LatticeAlternatingBinaryTraj::set_multiplier_for_bias_in_layer(int layer, double multiplier) {
        for (auto l: _lattices) {
            l.second->set_multiplier_for_bias_in_layer(layer, multiplier);
        };
    };
    
    // Add connections
    void LatticeAlternatingBinaryTraj::add_conn(int layer1, int x1, int layer2, int x2) {
        for (auto l: _lattices) {
            l.second->add_conn(layer1,x1,layer2,x2);
        };
    };
    void LatticeAlternatingBinaryTraj::add_conn(int layer1, int x1, int y1, int layer2, int x2, int y2) {
        for (auto l: _lattices) {
            l.second->add_conn(layer1,x1,y1,layer2,x2,y2);
        };
    };
    void LatticeAlternatingBinaryTraj::add_conn(int layer1, int x1, int y1, int z1, int layer2, int x2, int y2, int z2) {
        for (auto l: _lattices) {
            l.second->add_conn(layer1,x1,y1,z1,layer2,x2,y2,z2);
        };
    };
};
