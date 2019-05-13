#include "../include/dblz_bits/center.hpp"

// Other headers
#include "../include/dblz_bits/general.hpp"
#include "../include/dblz_bits/species.hpp"

#include <iostream>
#include <fstream>

/************************************
 * Namespace for dblz
 ************************************/

namespace dblz {
    
    /****************************************
     Center
     ****************************************/
    
    // ***************
    // MARK: - Constructor
    // ***************
    
    Center::Center(int layer, Sptr species, double val) {
        _val = val;
        _val_new = 0.0;
        _layer = layer;
        _species = species;
    };
    Center::Center(const Center& other) {
        _copy(other);
    };
    Center::Center(Center&& other) {
        _move(other);
    };
    Center& Center::operator=(const Center& other) {
        if (this != &other) {
            _clean_up();
            _copy(other);
        };
        return *this;
    };
    Center& Center::operator=(Center&& other) {
        if (this != &other) {
            _clean_up();
            _move(other);
        };
        return *this;
    };
    Center::~Center()
    {
        _clean_up();
    };
    void Center::_clean_up() {
        // Nothing....
    };
    void Center::_copy(const Center& other) {
        _val = other._val;
        _val_new = other._val_new;
        _layer = other._layer;
        _species = other._species;
    };
    void Center::_move(Center& other) {
        _val = other._val;
        _val_new = other._val_new;
        _layer = other._layer;
        _species = other._species;

        // Reset the other
        other._val = 0.5;
        other._val_new = 0.0;
        other._layer = 0;
        other._species = nullptr;
    };
    
    // ***************
    // MARK: - Getters/setters
    // ***************
    
    std::string Center::get_name() const {
        return _species->get_name() + "_layer_" + pad_str(_layer, 2);
    };
    
    void Center::set_val(double val) {
        _val = val;
    };
    double Center::get_val() const {
        return _val;
    };
    
    void Center::reset_val_new() {
        _val_new = 0.0;
    };
    void Center::set_val_new(double val_new) {
        _val_new = val_new;
    };
    void Center::increment_val_new(double inc) {
        _val_new += inc;
    };
    double Center::get_val_new() const {
        return _val_new;
    };
    
    void Center::slide(double sliding_factor) {
        // Slide
        _val = (1.0 - sliding_factor) * _val + sliding_factor * _val_new;
        
        // Reset
        _val_new = 0.0;
    };
    
    // ***************
    // MARK: - Get layer/species
    // ***************
    
    int Center::get_layer() const {
        return _layer;
    };
    Sptr Center::get_species() const {
        return _species;
    };
};
