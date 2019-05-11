#include "../include/dblz_bits/lattice_traj_alternating_binary.hpp"

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
	LatticeTrajAlternatingBinary
	****************************************/

	/********************
	Constructor
	********************/

    LatticeTrajAlternatingBinary::LatticeTrajAlternatingBinary(int no_dims, int box_length, std::vector<Sptr> species_visible) : LatticeTraj(no_dims, box_length, species_visible)
	{
        _lattices_alternating_binary[0] = std::make_shared<LatticeAlternatingBinary>(no_dims,box_length,species_visible);
        _lattices[0] = _lattices_alternating_binary[0];
        
        // Set no timesteps/timepoints
        set_no_timesteps(0,0);
	};
    LatticeTrajAlternatingBinary::LatticeTrajAlternatingBinary(const LatticeTrajAlternatingBinary& other) : LatticeTraj(other) {
		_copy(other);
	};
    LatticeTrajAlternatingBinary::LatticeTrajAlternatingBinary(LatticeTrajAlternatingBinary&& other) : LatticeTraj(std::move(other)) {
		_move(other);
	};
	LatticeTrajAlternatingBinary& LatticeTrajAlternatingBinary::operator=(const LatticeTrajAlternatingBinary& other) {
		if (this != &other) {
			_clean_up();
            LatticeTraj::operator=(other);
			_copy(other);
		};
		return *this;
	};
	LatticeTrajAlternatingBinary& LatticeTrajAlternatingBinary::operator=(LatticeTrajAlternatingBinary&& other) {
		if (this != &other) {
			_clean_up();
            LatticeTraj::operator=(std::move(other));
			_move(other);
		};
		return *this;
	};
	LatticeTrajAlternatingBinary::~LatticeTrajAlternatingBinary() {
		_clean_up();
	};

	void LatticeTrajAlternatingBinary::_clean_up() {
        // Nothing...
	};
	void LatticeTrajAlternatingBinary::_copy(const LatticeTrajAlternatingBinary& other) {
        for (auto l: other._lattices_alternating_binary) {
            _lattices_alternating_binary[l.first] = std::make_shared<LatticeAlternatingBinary>(*l.second);
            _lattices[l.first] = _lattices_alternating_binary.at(l.first);
        };
    };
	void LatticeTrajAlternatingBinary::_move(LatticeTrajAlternatingBinary& other) {
        _lattices = other._lattices;
        _lattices_alternating_binary = other._lattices_alternating_binary;

		// Reset other
        other._lattices.clear();
        other._lattices_alternating_binary.clear();
    };

    // ***************
    // MARK: - No timesteps
    // ***************
    
    void LatticeTrajAlternatingBinary::set_no_timesteps(int timepoint_start, int no_timesteps) {
        if (no_timesteps < 0) {
            std::cerr << ">>> LatticeTrajAlternatingBinary::set_no_timesteps <<< Error: no timepoints must be > 0!" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        // Save at least one lattice
        auto itl = _lattices_alternating_binary.end();
        itl--;
        auto lptr = itl->second;
        
        // Delete
        auto it = _lattices_alternating_binary.begin();
        while (it != _lattices_alternating_binary.end()) {
            if (it->first < timepoint_start || it->first > timepoint_start + no_timesteps) {
                // Delete
                it = _lattices_alternating_binary.erase(it);
            } else {
                it++;
            };
        };
        auto it2 = _lattices.begin();
        while (it2 != _lattices.end()) {
            if (it2->first < timepoint_start || it2->first > timepoint_start + no_timesteps) {
                // Delete
                it2 = _lattices.erase(it2);
            } else {
                it2++;
            };
        };
        
        // Add as needed
        for (auto timepoint=timepoint_start; timepoint<=timepoint_start + no_timesteps; timepoint++) {
            auto itf = _lattices_alternating_binary.find(timepoint);
            if (itf == _lattices_alternating_binary.end()) {
                // Does not exist; need to create
                _lattices_alternating_binary[timepoint] = std::make_shared<LatticeAlternatingBinary>(*lptr);
                _lattices[timepoint] = _lattices_alternating_binary.at(timepoint);
                
                // Set correct ixns
                for (auto pr1: _bias_dict) {
                    for (auto pr2: pr1.second) {
                        _lattices_alternating_binary[timepoint]->set_bias_of_layer(pr1.first, pr2.first, pr2.second->get_ixn_param_at_timepoint(timepoint));
                    };
                };
                for (auto pr1: _o2_ixn_dict) {
                    for (auto pr2: pr1.second) {
                        for (auto pr3: pr2.second) {
                            for (auto pr4: pr3.second) {
                                _lattices_alternating_binary[timepoint]->set_ixn_between_layers(pr1.first, pr2.first, pr3.first, pr4.first, pr4.second->get_ixn_param_at_timepoint(timepoint));
                            };
                        };
                    };
                };
            };
        };
    };
    
    // ***************
    // MARK: - Get lattice
    // ***************
    
    std::shared_ptr<LatticeAlternatingBinary> LatticeTrajAlternatingBinary::get_lattice_alternating_binary_at_timepoint(int timepoint) const {
        return _lattices_alternating_binary.at(timepoint);
    };
};
