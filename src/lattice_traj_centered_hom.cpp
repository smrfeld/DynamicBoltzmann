#include "../include/dblz_bits/lattice_traj_centered_hom.hpp"

// Other headers
#include "../include/dblz_bits/general.hpp"
#include "../include/dblz_bits/species.hpp"
#include "../include/dblz_bits/moment_diff.hpp"
#include "../include/dblz_bits/ixn_param.hpp"
#include "../include/dblz_bits/ixn_param_traj.hpp"
#include "../include/dblz_bits/lattice_centered_hom.hpp"
#include "../include/dblz_bits/adjoint.hpp"
#include "../include/dblz_bits/center_traj.hpp"
#include "../include/dblz_bits/center.hpp"
#include "../include/dblz_bits/fwds/fwds_center_traj.hpp"

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
	LatticeTrajCenteredHom
	****************************************/

	/********************
	Constructor
	********************/

    LatticeTrajCenteredHom::LatticeTrajCenteredHom(int no_dims, int box_length, std::vector<Sptr> species_visible, std::vector<CTptr> center_trajs) : LatticeTraj(no_dims, box_length, species_visible, false)
	{
        _center_trajs[0] = center_trajs;
        
        // Form centers at timepoint 0
        std::vector<Cptr> centers;
        for (auto center_traj: _center_trajs.at(0)) {
            centers.push_back(center_traj->get_center_at_timepoint(0));
        };
        
        // Make initial lattice
        _lattices_centered_hom[0] = std::make_shared<LatticeCenteredHom>(no_dims,box_length,species_visible,centers);
        _lattices[0] = _lattices_centered_hom[0];
        
        // Set no timesteps/timepoints
        set_no_timesteps(0,0);
	};
    LatticeTrajCenteredHom::LatticeTrajCenteredHom(const LatticeTrajCenteredHom& other) : LatticeTraj(other) {
		_copy(other);
	};
    LatticeTrajCenteredHom::LatticeTrajCenteredHom(LatticeTrajCenteredHom&& other) : LatticeTraj(std::move(other)) {
		_move(other);
	};
	LatticeTrajCenteredHom& LatticeTrajCenteredHom::operator=(const LatticeTrajCenteredHom& other) {
		if (this != &other) {
			_clean_up();
            LatticeTraj::operator=(other);
			_copy(other);
		};
		return *this;
	};
	LatticeTrajCenteredHom& LatticeTrajCenteredHom::operator=(LatticeTrajCenteredHom&& other) {
		if (this != &other) {
			_clean_up();
            LatticeTraj::operator=(std::move(other));
			_move(other);
		};
		return *this;
	};
	LatticeTrajCenteredHom::~LatticeTrajCenteredHom() {
		_clean_up();
	};

	void LatticeTrajCenteredHom::_clean_up() {
        // Nothing...
	};
	void LatticeTrajCenteredHom::_copy(const LatticeTrajCenteredHom& other) {
        for (auto l: other._lattices_centered_hom) {
            _lattices_centered_hom[l.first] = std::make_shared<LatticeCenteredHom>(*l.second);
            _lattices[l.first] = _lattices_centered_hom.at(l.first);
        };
        _center_trajs = other._center_trajs;
    };
	void LatticeTrajCenteredHom::_move(LatticeTrajCenteredHom& other) {
        _lattices = other._lattices;
        _lattices_centered_hom = other._lattices_centered_hom;
        _center_trajs = other._center_trajs;

		// Reset other
        other._lattices.clear();
        other._lattices_centered_hom.clear();
        other._center_trajs.clear();
    };

    // ***************
    // MARK: - No timesteps
    // ***************
    
    void LatticeTrajCenteredHom::set_no_timesteps(int timepoint_start, int no_timesteps) {
        if (no_timesteps < 0) {
            std::cerr << ">>> LatticeTrajCenteredHom::set_no_timesteps <<< Error: no timepoints must be > 0!" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        // Save at least one lattice
        auto itl = _lattices_centered_hom.end();
        itl--;
        auto lptr = itl->second;
        
        // Delete
        auto it = _lattices_centered_hom.begin();
        while (it != _lattices_centered_hom.end()) {
            if (it->first < timepoint_start || it->first > timepoint_start + no_timesteps) {
                // Delete
                it = _lattices_centered_hom.erase(it);
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
            auto itf = _lattices_centered_hom.find(timepoint);
            if (itf == _lattices_centered_hom.end()) {
                // Does not exist; need to create
                
                // Create lattice
                _lattices_centered_hom[timepoint] = std::make_shared<LatticeCenteredHom>(*lptr);
                _lattices[timepoint] = _lattices_centered_hom.at(timepoint);
                
                // Clear existing centers
                _lattices_centered_hom.at(timepoint)->clear_all_centers();
                
                // Set correct centers
                for (auto pr1: _center_trajs) {
                    std::vector<Cptr> centers;
                    for (auto center_traj: pr1.second) {
                        centers.push_back(center_traj->get_center_at_timepoint(timepoint));
                    };
                   
                    _lattices_centered_hom.at(timepoint)->set_centers_in_layer(pr1.first, centers);
                };
                
                // Clear existing ixns
                _lattices_centered_hom.at(timepoint)->clear_all_biases_and_ixns();
                
                // Set correct ixns
                for (auto pr1: _bias_dict) {
                    for (auto pr2: pr1.second) {
                        _lattices_centered_hom[timepoint]->set_bias_of_layer(pr1.first, pr2.first, pr2.second->get_ixn_param_at_timepoint(timepoint));
                    };
                };
                for (auto pr1: _o2_ixn_dict) {
                    for (auto pr2: pr1.second) {
                        for (auto pr3: pr2.second) {
                            for (auto pr4: pr3.second) {
                                _lattices_centered_hom[timepoint]->set_ixn_between_layers(pr1.first, pr2.first, pr3.first, pr4.first, pr4.second->get_ixn_param_at_timepoint(timepoint));
                            };
                        };
                    };
                };
            };
        };
    };
    
    // Add layer
    void LatticeTrajCenteredHom::add_layer(int layer, int box_length, std::vector<Sptr> species, std::vector<CTptr> center_trajs) {
        _center_trajs[layer] = center_trajs;
        
        // Go through lattices at all times
        for (auto l: _lattices_centered_hom) {
            
            // Form centers at timepoint
            std::vector<Cptr> centers;
            for (auto center_traj: center_trajs) {
                centers.push_back(center_traj->get_center_at_timepoint(l.first));
            };
            
            // Add layer
            l.second->add_layer(layer, box_length, species, centers);
        };
    };
    void LatticeTrajCenteredHom::add_layer(int layer, int box_length, std::vector<Sptr> species) {
        std::cerr << ">>> LatticeTrajCenteredHom::add_layer <<< must provide center trajs when adding layer" << std::endl;
        exit(EXIT_FAILURE);
    };

    // ***************
    // MARK: - Get lattice
    // ***************
    
    std::shared_ptr<LatticeCenteredHom> LatticeTrajCenteredHom::get_lattice_centered_hom_at_timepoint(int timepoint) const {
        return _lattices_centered_hom.at(timepoint);
    };
};
