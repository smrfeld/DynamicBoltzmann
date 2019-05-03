#include "../include/dblz_bits/center_traj.hpp"

// Other headers
#include "../include/dblz_bits/general.hpp"
#include "../include/dblz_bits/species.hpp"
#include "../include/dblz_bits/center.hpp"

#include <iostream>
#include <fstream>

/************************************
 * Namespace for dblz
 ************************************/

namespace dblz {
    
    /****************************************
     CenterTraj
     ****************************************/
    
    // ***************
    // MARK: - Constructor
    // ***************
    
    CenterTraj::CenterTraj(int layer, Sptr species) {
        _layer = layer;
        _species = species;
        
        set_no_timesteps(0);
    };
    CenterTraj::CenterTraj(const CenterTraj& other) {
        _copy(other);
    };
    CenterTraj::CenterTraj(CenterTraj&& other) {
        _move(other);
    };
    CenterTraj& CenterTraj::operator=(const CenterTraj& other) {
        if (this != &other) {
            _clean_up();
            _copy(other);
        };
        return *this;
    };
    CenterTraj& CenterTraj::operator=(CenterTraj&& other) {
        if (this != &other) {
            _clean_up();
            _move(other);
        };
        return *this;
    };
    CenterTraj::~CenterTraj()
    {
        _clean_up();
    };
    void CenterTraj::_clean_up() {
        // Nothing....
    };
    void CenterTraj::_copy(const CenterTraj& other) {
        _no_timesteps = other._no_timesteps;
        _no_timepoints = other._no_timepoints;
        _centers = other._centers;
        _layer = other._layer;
        _species = other._species;
    };
    void CenterTraj::_move(CenterTraj& other) {
        _no_timesteps = other._no_timesteps;
        _no_timepoints = other._no_timepoints;
        _centers = other._centers;
        _layer = other._layer;
        _species = other._species;

        // Reset the other
        other._no_timesteps = 0;
        other._no_timepoints = 0;
        other._centers.clear();
        other._layer = 0;
        other._species = nullptr;
    };
    
    // ***************
    // MARK: - Timesteps
    // ***************
    
    int CenterTraj::get_no_timesteps() const {
        return _no_timesteps;
    };
    void CenterTraj::set_no_timesteps(int no_timesteps) {
        _no_timesteps = no_timesteps;
        _no_timepoints = _no_timesteps + 1;
        
        while (_centers.size() < _no_timepoints) {
            _centers.push_back(std::make_shared<Center>(_layer,_species));
        };
        while (_centers.size() > _no_timepoints) {
            _centers.pop_back();
        };
    };
    
    // ***************
    // MARK: - Getters/setters
    // ***************
    
    std::string CenterTraj::get_name() const {
        return _species->get_name() + "_layer_" + pad_str(_layer, 2);
    };
    
    Cptr CenterTraj::get_center_at_timepoint(int timepoint) const {
        return _centers.at(timepoint);
    };

    double CenterTraj::get_val_at_timepoint(int timepoint) const {
        return _centers.at(timepoint)->get_val();
    };
    
    // Time derivative
    double CenterTraj::get_deriv_at_timepoint(int timepoint, double dt) const {
        return (_centers.at(timepoint+1)->get_val() - _centers.at(timepoint)->get_val()) / dt;
    };
    
    // ***************
    // MARK: - Get layer/species
    // ***************
    
    int CenterTraj::get_layer() const {
        return _layer;
    };
    Sptr CenterTraj::get_species() const {
        return _species;
    };
};
