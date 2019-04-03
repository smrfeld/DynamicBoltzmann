#include "../include/dblz_bits/lattice_traj_1d_fully_visible.hpp"

// Other headers
#include "../include/dblz_bits/general.hpp"
#include "../include/dblz_bits/species.hpp"
#include "../include/dblz_bits/moment.hpp"
#include "../include/dblz_bits/ixn_param.hpp"
#include "../include/dblz_bits/ixn_param_traj.hpp"
#include "../include/dblz_bits/lattice_1d_fully_visible.hpp"
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
     LatticeTraj1DFullyVisible
     ****************************************/
    
    /********************
     Constructor
     ********************/
    
    LatticeTraj1DFullyVisible::LatticeTraj1DFullyVisible(int box_length, std::vector<Sptr> species_visible, bool conn_2, bool conn_3)
    {
        _conn_2 = conn_2;
        _conn_3 = conn_3;
        _lattices[0] = std::make_shared<Lattice1DFullyVisible>(box_length,species_visible,conn_2,conn_3);
        
        // Set no timesteps/timepoints
        set_no_timesteps(0,0);
    };
    LatticeTraj1DFullyVisible::LatticeTraj1DFullyVisible(const LatticeTraj1DFullyVisible& other) {
        _copy(other);
    };
    LatticeTraj1DFullyVisible::LatticeTraj1DFullyVisible(LatticeTraj1DFullyVisible&& other) {
        _move(other);
    };
    LatticeTraj1DFullyVisible& LatticeTraj1DFullyVisible::operator=(const LatticeTraj1DFullyVisible& other) {
        if (this != &other) {
            _clean_up();
            _copy(other);
        };
        return *this;
    };
    LatticeTraj1DFullyVisible& LatticeTraj1DFullyVisible::operator=(LatticeTraj1DFullyVisible&& other) {
        if (this != &other) {
            _clean_up();
            _move(other);
        };
        return *this;
    };
    LatticeTraj1DFullyVisible::~LatticeTraj1DFullyVisible() {
        _clean_up();
    };
    
    void LatticeTraj1DFullyVisible::_clean_up() {
        // Nothing...
    };
    void LatticeTraj1DFullyVisible::_copy(const LatticeTraj1DFullyVisible& other) {
        for (auto l: other._lattices) {
            _lattices[l.first] = std::make_shared<Lattice1DFullyVisible>(*l.second);
        };
        _conn_2 = other._conn_2;
        _conn_3 = other._conn_3;
        _ixn_param_trajs = other._ixn_param_trajs;
        _bias_dict = other._bias_dict;
        _o2_ixn_dict = other._o2_ixn_dict;
    };
    void LatticeTraj1DFullyVisible::_move(LatticeTraj1DFullyVisible& other) {
        _lattices = other._lattices;
        _conn_2 = other._conn_2;
        _conn_3 = other._conn_3;
        _ixn_param_trajs = other._ixn_param_trajs;
        _bias_dict = other._bias_dict;
        _o2_ixn_dict = other._o2_ixn_dict;
        
        // Reset other
        other._lattices.clear();
        other._conn_2 = false;
        other._conn_3 = false;
        other._ixn_param_trajs.clear();
        other._bias_dict.clear();
        other._o2_ixn_dict.clear();
    };
    
    // ***************
    // MARK: - [Private] Add ixn param traj
    // ***************
    
    void LatticeTraj1DFullyVisible::_add_ixn_param_traj(ITptr ixn_param_traj) {
        auto it = std::find(_ixn_param_trajs.begin(),_ixn_param_trajs.end(),ixn_param_traj);
        if (it == _ixn_param_trajs.end()) {
            _ixn_param_trajs.push_back(ixn_param_traj);
        };
    };
    
    // ***************
    // MARK: - No timesteps
    // ***************
    
    void LatticeTraj1DFullyVisible::set_no_timesteps(int timepoint_start, int no_timesteps) {
        if (no_timesteps < 0) {
            std::cerr << ">>> LatticeTraj1DFullyVisible::set_no_timesteps <<< Error: no timepoints must be > 0!" << std::endl;
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
                _lattices[timepoint] = std::make_shared<Lattice1DFullyVisible>(*lptr);
                
                // Set correct ixns
                for (auto pr: _bias_dict) {
                    _lattices[timepoint]->set_bias(pr.first, pr.second->get_ixn_param_at_timepoint(timepoint));
                };
                for (auto pr1: _o2_ixn_dict) {
                    for (auto pr2: pr1.second) {
                        _lattices[timepoint]->set_ixn_2(pr1.first, pr2.first, pr2.second->get_ixn_param_at_timepoint(timepoint));
                    };
                };
                for (auto pr1: _o3_ixn_dict) {
                    for (auto pr2: pr1.second) {
                        for (auto pr3: pr2.second) {
                            _lattices[timepoint]->set_ixn_3(pr1.first, pr2.first, pr3.first, pr3.second->get_ixn_param_at_timepoint(timepoint));
                        };
                    };
                };
            };
        };
    };
    
    // ***************
    // MARK: - Getters
    // ***************
    
    int LatticeTraj1DFullyVisible::get_box_length() const {
        return _lattices.begin()->second->get_box_length();
    };
    
    // Get all ixns
    const std::vector<ITptr>& LatticeTraj1DFullyVisible::get_all_ixn_param_trajs() const {
        return _ixn_param_trajs;
    };
    
    // ***************
    // MARK: - Markov chains
    // ***************
    
    int LatticeTraj1DFullyVisible::get_no_markov_chains(MCType type) const {
        return _lattices.begin()->second->get_no_markov_chains(type);
    };
    
    void LatticeTraj1DFullyVisible::set_no_markov_chains(MCType type, int no_markov_chains) {
        for (auto l: _lattices) {
            l.second->set_no_markov_chains(type, no_markov_chains);
        };
    };
    
    // ***************
    // MARK: - Get lattice
    // ***************
    
    std::shared_ptr<Lattice1DFullyVisible> LatticeTraj1DFullyVisible::get_lattice_at_timepoint(int timepoint) const {
        return _lattices.at(timepoint);
    };
    
    // ***************
    // MARK: - Setup all timepoints
    // ***************
    
    // Biases
    void LatticeTraj1DFullyVisible::set_bias(Sptr sp, ITptr bias) {
        _bias_dict[sp] = bias;
        
        _add_ixn_param_traj(bias);
        
        for (auto l: _lattices) {
            l.second->set_bias(sp, bias->get_ixn_param_at_timepoint(l.first));
        };
    };
    
    // Ixns
    void LatticeTraj1DFullyVisible::set_ixn_2(Sptr sp1, Sptr sp2, ITptr ixn) {
        _o2_ixn_dict[sp1][sp2] = ixn;
        _o2_ixn_dict[sp2][sp1] = ixn;
        
        _add_ixn_param_traj(ixn);
        
        for (auto l: _lattices) {
            l.second->set_ixn_2(sp1, sp2, ixn->get_ixn_param_at_timepoint(l.first));
        };
    };
    void LatticeTraj1DFullyVisible::set_ixn_3(Sptr sp1, Sptr sp2, Sptr sp3, ITptr ixn) {
        _o3_ixn_dict[sp1][sp2][sp3] = ixn;
        _o3_ixn_dict[sp3][sp2][sp1] = ixn;
        
        _add_ixn_param_traj(ixn);
        
        for (auto l: _lattices) {
            l.second->set_ixn_3(sp1, sp2, sp3, ixn->get_ixn_param_at_timepoint(l.first));
        };
    };

};
