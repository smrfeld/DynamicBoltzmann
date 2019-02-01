#include "../include/bmla_bits/lattice.hpp"

// Other headers
#include "../include/bmla_bits/general.hpp"
#include "../include/bmla_bits/species.hpp"
#include "../include/bmla_bits/moment.hpp"
#include "../include/bmla_bits/ixn_param.hpp"

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

namespace bmla {

	/****************************************
	Lattice
	****************************************/

	/********************
	Constructor
	********************/

	Lattice::Lattice(int no_dims, int box_length, std::vector<Sptr> species_visible, bool batch_norm_mode)
	{
		if (no_dims != 1 && no_dims != 2 && no_dims != 3) {
			std::cerr << "ERROR: only dimensions 1,2,3 are supported for Lattice." << std::endl;
			exit(EXIT_FAILURE);
		};
        
		_no_dims = no_dims;
		_box_length = box_length;
        _no_layers = 0;
        
        // Batch norm
        _bn_mode = batch_norm_mode;
        _bn_eps = 1.0e-8;

        // Set no chains
        set_no_markov_chains(MCType::AWAKE, 1);
        set_no_markov_chains(MCType::ASLEEP, 1);
    
		// Visible layer
        add_layer(0, _box_length, species_visible, nullptr, nullptr);
	};
	Lattice::Lattice(const Lattice& other) {
		_copy(other);
	};
	Lattice::Lattice(Lattice&& other) {
		_move(other);
	};
	Lattice& Lattice::operator=(const Lattice& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	Lattice& Lattice::operator=(Lattice&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
	};
	Lattice::~Lattice() {
		_clean_up();
	};

	void Lattice::_clean_up() {
        // Nothing...
	};
	void Lattice::_copy(const Lattice& other) {
		_no_dims = other._no_dims;
		_box_length = other._box_length;
        _no_markov_chains = other._no_markov_chains;
        _no_layers = other._no_layers;
        
        _mc_chains = other._mc_chains;
        _mc_chains_act = other._mc_chains_act;
        _no_units_per_layer = other._no_units_per_layer;
        
        _lookup_1 = other._lookup_1;
        _lookup_2 = other._lookup_2;
        _lookup_3 = other._lookup_3;
        _rlookup = other._rlookup;
        
        _species_possible_map = other._species_possible_map;
        _species_possible_vec = other._species_possible_vec;
        
        _adj = other._adj;
        
        _all_ixns = other._all_ixns;
        _bias_dict = other._bias_dict;
        _o2_ixn_dict = other._o2_ixn_dict;
        _o2_mults = other._o2_mults;
      
        _bn_mode = other._bn_mode;
        _bn_beta = other._bn_beta;
        _bn_gamma = other._bn_gamma;
        _bn_beta_bar = other._bn_beta_bar;
        _bn_gamma_bar = other._bn_gamma_bar;
        _bn_means = other._bn_means;
        _bn_vars = other._bn_vars;
        _bn_eps = other._bn_eps;
    };
	void Lattice::_move(Lattice& other) {
        _no_dims = other._no_dims;
        _box_length = other._box_length;
        _no_markov_chains = other._no_markov_chains;
        _no_layers = other._no_layers;

        _mc_chains = other._mc_chains;
        _mc_chains_act = other._mc_chains_act;
        _no_units_per_layer = other._no_units_per_layer;
        
        _lookup_1 = other._lookup_1;
        _lookup_2 = other._lookup_2;
        _lookup_3 = other._lookup_3;
        _rlookup = other._rlookup;
        
        _species_possible_map = other._species_possible_map;
        _species_possible_vec = other._species_possible_vec;
        
        _adj = other._adj;
        
        _all_ixns = other._all_ixns;
        _bias_dict = other._bias_dict;
        _o2_ixn_dict = other._o2_ixn_dict;
        _o2_mults = other._o2_mults;
        
        _bn_mode = other._bn_mode;
        _bn_beta = other._bn_beta;
        _bn_gamma = other._bn_gamma;
        _bn_beta_bar = other._bn_beta_bar;
        _bn_gamma_bar = other._bn_gamma_bar;
        _bn_means = other._bn_means;
        _bn_vars = other._bn_vars;
        _bn_eps = other._bn_eps;

		// Reset other
		other._no_dims = 0;
		other._box_length = 0;
        other._no_markov_chains.clear();
        other._no_layers = 0;
        
        other._mc_chains.clear();
        other._mc_chains_act.clear();
        other._no_units_per_layer.clear();
        
		other._lookup_1.clear();
		other._lookup_2.clear();
		other._lookup_3.clear();
        other._rlookup.clear();
        
        other._species_possible_map.clear();
        other._species_possible_vec.clear();
        
        other._adj.clear();
        
        other._all_ixns.clear();
        other._bias_dict.clear();
        other._o2_ixn_dict.clear();
        other._o2_mults.clear();
        
        other._bn_mode = true;
        other._bn_beta.clear();
        other._bn_gamma.clear();
        other._bn_beta_bar.clear();
        other._bn_gamma_bar.clear();
        other._bn_means.clear();
        other._bn_vars.clear();
        other._bn_eps = 0.0;
	};
    
    /****************************************
    PRIVATE METHODS
     ****************************************/

    // Lookup a site iterator from x,y,z
    int Lattice::_look_up_unit(int layer, int x) const {
        auto it = _lookup_1.find(layer);
        if (it != _lookup_1.end()) {
            auto it2 = it->second.find(x);
            if (it2 != it->second.end()) {
                return it2->second;
            };
        };
        
        std::cerr << ">>> Lattice::_look_up_unit <<< could not find layer: " << layer << " x: " << x << std::endl;
        exit(EXIT_FAILURE);
    };
    int Lattice::_look_up_unit(int layer, int x, int y) const {
        auto it = _lookup_2.find(layer);
        if (it != _lookup_2.end()) {
            auto it2 = it->second.find(x);
            if (it2 != it->second.end()) {
                auto it3 = it2->second.find(y);
                if (it3 != it2->second.end()) {
                    return it3->second;
                };
            };
        };
        
        std::cerr << ">>> Lattice::_look_up_unit <<< could not find layer: " << layer << " x: " << x << " y: " << y << std::endl;
        exit(EXIT_FAILURE);
    };
    int Lattice::_look_up_unit(int layer, int x, int y, int z) const {
        auto it = _lookup_3.find(layer);
        if (it != _lookup_3.end()) {
            auto it2 = it->second.find(x);
            if (it2 != it->second.end()) {
                auto it3 = it2->second.find(y);
                if (it3 != it2->second.end()) {
                    auto it4 = it3->second.find(z);
                    if (it4 != it3->second.end()) {
                        return it4->second;
                    };
                };
            };
        };
        
        std::cerr << ">>> Lattice::_look_up_unit <<< could not find layer: " << layer << " x: " << x << " y: " << y << " z: " << z << std::endl;
        exit(EXIT_FAILURE);
    };
    
    std::vector<int> Lattice::_look_up_pos(int layer, int idx) const {
        auto it = _rlookup.find(layer);
        if (it != _rlookup.end()) {
            auto it2 = it->second.find(idx);
            if (it2 != it->second.end()) {
                return it2->second;
            };
        };

        return std::vector<int>();
    };
    
    /****************************************
    PUBLIC METHODS
     ****************************************/

    /********************
     Getters
     ********************/
    
    int Lattice::get_no_dims() const {
        return _no_dims;
    };
    int Lattice::get_box_length() const {
        return _box_length;
    };
    
    int Lattice::get_no_units_in_layer(int layer) const {
        return _no_units_per_layer.at(layer);
    };
    
    int Lattice::get_no_layers() const {
        return _no_layers;
    };
    
    /********************
     Markov chains
     ********************/
    
    int Lattice::get_no_markov_chains(MCType type) const {
        return _no_markov_chains.at(type);
    };

    void Lattice::set_no_markov_chains(MCType type, int no_markov_chains) {
        _no_markov_chains[type] = no_markov_chains;
        
        // Add chains
        int no_units;
        if (_mc_chains[type].size() < _no_markov_chains[type]) {
            for (auto i_chain=_mc_chains[type].size(); i_chain < _no_markov_chains[type]; i_chain++) {
                for (auto layer=0; layer<_no_layers; layer++) {
                    no_units = get_no_units_in_layer(layer);
                    for (auto sp: _species_possible_vec[layer]) {
                        _mc_chains[type][i_chain][layer][sp] = arma::vec(no_units,arma::fill::randu);
                        _mc_chains_act[type][i_chain][layer][sp] = arma::vec(no_units,arma::fill::randu);
                    };
                };
            };
        };
        
        // Remove chains
        if (_mc_chains[type].size() > _no_markov_chains[type]) {
            for (auto i_chain=_mc_chains[type].size(); i_chain > _no_markov_chains[type]; i_chain--) {
                _mc_chains[type].erase(i_chain);
                _mc_chains_act[type].erase(i_chain);
            };
        };
    };
    
    /********************
     Add a layer
     ********************/
    
    void Lattice::add_layer(int layer, int box_length, std::vector<Sptr> species) {
        if (_bn_mode) {
            std::cerr << ">>> Lattice::add_layer <<< in batch norm mode; must specify scale and shift parameters" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        add_layer(layer,box_length,species,nullptr,nullptr);
    };
    
    void Lattice::add_layer(int layer, int box_length, std::vector<Sptr> species, Iptr beta, Iptr gamma) {
        if (!_bn_mode && beta != nullptr && gamma != nullptr) {
            std::cerr << ">>> Lattice::add_layer <<< NOT in batch norm mode; cannot specify scale and shift parameters" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        if (layer != _no_layers) {
            std::cerr << ">>> Lattice::add_layer <<< error: next layer must be: " << _no_layers << " not: " << layer << std::endl;
            exit(EXIT_FAILURE);
        };
        
        // Increment no layers
        _no_layers += 1;
        
        // Add random
        int no_units = pow(box_length,_no_dims);
        for (auto &chain: _no_markov_chains) {
            for (auto i_chain=0; i_chain<chain.second; i_chain++) {
                for (auto sp: species) {
                    _mc_chains[chain.first][i_chain][layer][sp] = arma::vec(no_units,arma::fill::randu);
                    _mc_chains_act[chain.first][i_chain][layer][sp] = arma::vec(no_units,arma::fill::randu);
                };
            };
        };
        
        // No units
        _no_units_per_layer[layer] = no_units;

        // Lookups
        if (_no_dims == 1) {
            for (auto x=1; x<=_box_length; x++) {
                _lookup_1[layer][x] = x-1;
                _rlookup[layer][x-1] = std::vector<int>({x});
            };
        } else if (_no_dims == 2) {
            int ctr=0;
            for (auto x=1; x<=_box_length; x++) {
                for (auto y=1; y<=_box_length; y++) {
                    _lookup_2[layer][x][y] = ctr;
                    _rlookup[layer][ctr++] = std::vector<int>({x,y});
                };
            };
        } else if (_no_dims == 3) {
            int ctr=0;
            for (auto x=1; x<=_box_length; x++) {
                for (auto y=1; y<=_box_length; y++) {
                    for (auto z=1; z<=_box_length; z++) {
                        _lookup_3[layer][x][y][z] = ctr;
                        _rlookup[layer][ctr++] = std::vector<int>({x,y,z});
                    };
                };
            };
        };
        
        // Add species possible
        for (auto sp: species) {
            _species_possible_map[layer][sp->get_name()] = sp;
            _species_possible_vec[layer].push_back(sp);
        };
        
        // Add adjacency matrix
        if (layer != 0) {
            int size_below = get_no_units_in_layer(layer-1);
            _adj[layer-1] = arma::mat(size_below,no_units,arma::fill::zeros);
        };
        
        // Batch norm
        if (_bn_mode) {
            if (layer != 0) {
              
                for (auto sp: _species_possible_vec[layer]) {
                    // Set params
                    _bn_beta[layer][sp] = beta;
                    _bn_gamma[layer][sp] = gamma;
                    
                    // Add to all
                    auto it1 = std::find(_all_ixns.begin(), _all_ixns.end(), beta);
                    if (it1 == _all_ixns.end()) {
                        _all_ixns.push_back(beta);
                    };
                    auto it2 = std::find(_all_ixns.begin(), _all_ixns.end(), gamma);
                    if (it2 == _all_ixns.end()) {
                        _all_ixns.push_back(gamma);
                    };

                    _bn_beta_bar[MCType::AWAKE][layer][sp] = arma::vec(no_units,arma::fill::zeros);
                    _bn_beta_bar[MCType::ASLEEP][layer][sp] = arma::vec(no_units,arma::fill::zeros);
                    _bn_gamma_bar[MCType::AWAKE][layer][sp] = arma::vec(no_units,arma::fill::zeros);
                    _bn_gamma_bar[MCType::ASLEEP][layer][sp] = arma::vec(no_units,arma::fill::zeros);

                    _bn_means[MCType::AWAKE][layer][sp] = arma::vec(no_units,arma::fill::zeros);
                    _bn_means[MCType::ASLEEP][layer][sp] = arma::vec(no_units,arma::fill::zeros);
                    _bn_vars[MCType::AWAKE][layer][sp] = arma::vec(no_units,arma::fill::zeros);
                    _bn_vars[MCType::ASLEEP][layer][sp] = arma::vec(no_units,arma::fill::zeros);
                };
            };
        };
    };

	/********************
	Helpers to setup all sites - Biases
	********************/

    // Biases
	void Lattice::add_bias_all_layers(Sptr sp, Iptr bias) {
        for (auto layer=0; layer<_no_layers; layer++) {
            add_bias_to_layer(layer, sp, bias);
        };
    };

    void Lattice::add_bias_to_layer(int layer, Sptr sp, Iptr bias) {
        _bias_dict[layer][sp].push_back(bias);
        
        // Add to all
        auto it = std::find(_all_ixns.begin(), _all_ixns.end(), bias);
        if (it == _all_ixns.end()) {
            _all_ixns.push_back(bias);
        };
    };

    // Ixns
    void Lattice::add_ixn_between_layers(int layer1, Sptr sp1, int layer2, Sptr sp2, Iptr ixn) {
        if (layer2 != layer1+1 && layer2 != layer1-1) {
            std::cerr << ">>> Lattice::add_ixn_between_layers <<< layer2 != layer1 +- 1; instead layer_2 = " << layer2 << " and layer1 = " << layer1 << std::endl;
            exit(EXIT_FAILURE);
        };
        
        // Add both ways
        _o2_ixn_dict[layer1][sp1][layer2][sp2].push_back(ixn);
        _o2_ixn_dict[layer2][sp2][layer1][sp1].push_back(ixn);

        // Add to all
        auto it = std::find(_all_ixns.begin(), _all_ixns.end(), ixn);
        if (it == _all_ixns.end()) {
            _all_ixns.push_back(ixn);
        };
    };
    
    // Set multiplier
    void Lattice::set_multiplier_between_layers(int from_layer, int to_layer, double multiplier) {
        if (to_layer != from_layer+1 && to_layer != from_layer-1) {
            std::cerr << ">>> Lattice::set_multiplier <<< to_layer != from_layer +- 1; instead layer_2 = " << to_layer << " and from_layer = " << from_layer << std::endl;
            exit(EXIT_FAILURE);
        };
        
        _o2_mults[from_layer][to_layer] = multiplier;
    };

    // Get ixns
    double Lattice::get_bias_in_layer(int layer, Sptr sp) const {
        auto it = _bias_dict.find(layer);
        double val = 0.0;
        if (it != _bias_dict.end()) {
            auto it2 = it->second.find(sp);
            if (it2 != it->second.end()) {
                for (auto ixn: it2->second) {
                    val += ixn->get_val();
                };
            };
        };
        return val;
    };
    double Lattice::get_ixn_between_layers(int from_layer, Sptr from_sp, int to_layer, Sptr to_sp) const {
        if (to_layer != from_layer+1 && to_layer != from_layer-1) {
            std::cerr << ">>> Lattice::add_ixn_between_layers <<< to_layer != from_layer +- 1; instead to_layer = " << to_layer << " and from_layer = " << from_layer << std::endl;
            exit(EXIT_FAILURE);
        };
        
        // Multiplier
        double mult = 1.0;
        auto itm = _o2_mults.find(from_layer);
        if (itm != _o2_mults.end()) {
            auto itm2 = itm->second.find(to_layer);
            if (itm2 != itm->second.end()) {
                mult = itm2->second;
                // std::cout << "Lattice::get_ixn_between_layers: from_layer = " << from_layer << " to to_layer = " << to_layer << " mult = " << mult << std::endl;
            };
        };
        
        double val = 0.0;
        auto it1 = _o2_ixn_dict.find(from_layer);
        if (it1 != _o2_ixn_dict.end()) {
            auto it2 = it1->second.find(from_sp);
            if (it2 != it1->second.end()) {
                auto it3 = it2->second.find(to_layer);
                if (it3 != it2->second.end()) {
                    auto it4 = it3->second.find(to_sp);
                    if (it4 != it3->second.end()) {
                        for (auto ixn: it4->second) {
                            val += ixn->get_val();
                        };
                    };
                };
            };
        };
        
        return mult * val;
    };

    // Get all ixns
    const std::vector<Iptr>& Lattice::get_all_ixn_params() const {
        return _all_ixns;
    };
    
	/********************
	Helpers to setup all sites - Visible-Visible ixns
	********************/

    void Lattice::add_conn(int layer1, int x1, int layer2, int x2) {
        int idx1 = _look_up_unit(layer1, x1);
        int idx2 = _look_up_unit(layer2, x2);
        _adj[layer1](idx1,idx2) = 1.0;
    };
    void Lattice::add_conn(int layer1, int x1, int y1, int layer2, int x2, int y2) {
        int idx1 = _look_up_unit(layer1, x1, y1);
        int idx2 = _look_up_unit(layer2, x2, y2);
        _adj[layer1](idx1,idx2) = 1.0;
    };
    void Lattice::add_conn(int layer1, int x1, int y1, int z1, int layer2, int x2, int y2, int z2) {
        int idx1 = _look_up_unit(layer1, x1, y1, z1);
        int idx2 = _look_up_unit(layer2, x2, y2, z2);
        _adj[layer1](idx1,idx2) = 1.0;
    };
    
	/********************
	Apply funcs to all units
	********************/

	// Clear the lattice
	void Lattice::set_empty_all_units(MCType chain, int i_chain) {
        for (auto layer=0; layer<_no_layers; layer++) {
            set_empty_all_units_in_layer(chain, i_chain, layer);
        };
	};
    void Lattice::set_empty_all_units_in_layer(MCType chain, int i_chain, int layer) {
        for (auto sp: _species_possible_vec[layer]) {
            _mc_chains[chain][i_chain][layer][sp].fill(0.0);
        };
    };
    void Lattice::set_empty_all_hidden_units(MCType chain, int i_chain) {
        for (auto layer=1; layer<_no_layers; layer++) {
            set_empty_all_units_in_layer(chain, i_chain, layer);
        };
    };

    // Random
    void Lattice::set_random_all_units(MCType chain, int i_chain, bool binary) {
        for (auto layer=0; layer<_no_layers; layer++) {
            set_random_all_units_in_layer(chain, i_chain, layer,binary);
        };
    };
    void Lattice::set_random_all_units_in_layer(MCType chain, int i_chain, int layer, bool binary) {
        // Random probs
        for (auto sp: _species_possible_vec.at(layer)) {
            _mc_chains[chain][i_chain][layer][sp].fill(arma::fill::randu);
        };
        
        if (binary) {
            binarize_all_units_in_layer(chain, i_chain, layer);
        };
    };
    void Lattice::set_random_all_hidden_units(MCType chain, int i_chain, bool binary) {
        for (auto layer=1; layer<_no_layers; layer++) {
            set_random_all_units_in_layer(chain,i_chain,layer,binary);
        };
    };

    // Binarize
    void Lattice::binarize_all_units(MCType chain, int i_chain) {
        for (auto layer=0; layer<_no_layers; layer++) {
            binarize_all_units_in_layer(chain,i_chain,layer);
        };
    };
    void Lattice::_binarize_all_units_in_layer(MCType chain, int i_chain, int layer, bool act) {
        // Random vec
        int no_units = get_no_units_in_layer(layer);
        auto r = arma::vec(no_units,arma::fill::randu);
        
        // Helpers
        
        // Running for evaluating propensities
        auto running = arma::vec(no_units,arma::fill::zeros);
        
        // flip vector:
        // 1 initially = rand is above current propensity
        // 0 = rand is below current propensity
        auto flip = arma::vec(no_units,arma::fill::ones);
        auto new_flip = flip;

        // _latt or _latt_act?
        if (act) {
            
            // Convert to propensities
            for (auto sp: _species_possible_vec.at(layer)) {
                running += _mc_chains_act[chain][i_chain][layer][sp];
                _mc_chains_act[chain][i_chain][layer][sp] = running;
            };
            
            // Evaluate
            // Initially: we (flip) are at zero, rand value is somewhere above
            // Next: we (flip) are at some propensity, rand value may be somewhere above or below
            // Finally: we (flip) are at max propensity = 1.0, rand value is below
            for (auto sp: _species_possible_vec.at(layer)) {
                // New flip vector
                new_flip = 0.5*sign(r - _mc_chains_act[chain][i_chain][layer][sp]) + 0.5;
                // flip -> new_flip
                // 1 -> 1 => not this species => (flip - new_flip) = 0
                // 1 -> 0 => yes this species => (flip - new_flip) = 1
                // 0 -> 0 => already assigned a previous species  => (flip - new_flip) = 0
                _mc_chains_act[chain][i_chain][layer][sp] = flip - new_flip;
                
                // Old flip vector = new flip vector
                flip = new_flip;
            };
            // Finally: if still above the final propensity, it is the empty species! This means all are zero, so no further adjustment needed!
            
        } else {
            
            // As above but with _latt
            for (auto sp: _species_possible_vec.at(layer)) {
                running += _mc_chains[chain][i_chain][layer][sp];
                _mc_chains[chain][i_chain][layer][sp] = running;
            };
            
            for (auto sp: _species_possible_vec.at(layer)) {
                new_flip = 0.5*sign(r - _mc_chains[chain][i_chain][layer][sp]) + 0.5;
                _mc_chains[chain][i_chain][layer][sp] = flip - new_flip;
                flip = new_flip;
            };
        };
    };
    void Lattice::binarize_all_units_in_layer(MCType chain, int i_chain, int layer) {
        _binarize_all_units_in_layer(chain,i_chain,layer,false);
    };
    void Lattice::binarize_all_hidden_units(MCType chain, int i_chain) {
        for (auto layer=1; layer<_no_layers; layer++) {
            binarize_all_units_in_layer(chain,i_chain,layer);
        };
    };

	/********************
	Write/read Latt to a file
	********************/

	void Lattice::write_layer_to_file(MCType chain, int i_chain, int layer, std::string fname, bool binary) const
	{
		std::ofstream f;
		f.open (fname);
		if (!f.is_open()) { // make sure we found it
			std::cerr << ">>> Lattice::write_layer_to_file <<< Error: could not open file: " << fname << " for writing" << std::endl;
			exit(EXIT_FAILURE);
		};
        
        // Go through all units
        int no_units = get_no_units_in_layer(layer);
        auto slatt = _mc_chains.at(chain).at(i_chain).at(layer);
        auto sps = _species_possible_vec.at(layer);
        std::vector<int> pos;
        if (binary) {

            // Binary

            // Go through all units
            for (auto unit=0; unit<no_units; unit++) {
                
                // Go through species possible
                for (auto sp: sps) {
                    if (abs(slatt.at(sp)(unit) - 1.0) < 1.0e-5) {
                        // Write pos
                        pos = _look_up_pos(layer,unit);
                        for (auto const &x: pos) {
                            f << x << " ";
                        };
                        
                        // Write species
                        f << sp->get_name() << "\n";
                        break;
                    };
                };
            };
                
        } else {
                
            // Not binary
            
            // Go through all units
            for (auto unit=0; unit<no_units; unit++) {
                
                // Write pos
                pos = _look_up_pos(layer,unit);
                for (auto const &x: pos) {
                    f << x << " ";
                };

                // Go through species possible
                for (auto sp: sps) {
                    // Write species
                    f << sp->get_name() << " " << slatt.at(sp)(unit) << " ";
                };
                f << "\n";
            };
        };
        
		f.close();
	};

	void Lattice::read_layer_from_file(MCType chain, int i_chain, int layer, std::string fname, bool binary)
	{
		// Clear the layer
        set_empty_all_units_in_layer(chain, i_chain, layer);

		// Open
		std::ifstream f;
		f.open(fname);
		if (!f.is_open()) { // make sure we found it
			std::cerr << ">>> Error: Lattice::read_layer_from_file <<< could not find file: " << fname << std::endl;
			exit(EXIT_FAILURE);
		};

        std::map<std::string, std::vector<int>> occs_to_write_binary;
        std::map<std::string, std::vector<std::pair<int,double>>> occs_to_write_prob;

        std::string x="",y="",z="";
        std::string sp="";
        std::string line;
        std::istringstream iss;
        int s;
        std::string prob="";
        double prob_val;

        if (_no_dims == 1 && binary) {
            while (getline(f,line)) {
                if (line == "") { continue; };
                iss = std::istringstream(line);
                iss >> x >> sp;
                s = _look_up_unit(layer,atoi(x.c_str()));
                occs_to_write_binary[sp].push_back(s);
                // Reset
                sp=""; x="";
            };
        } else if (_no_dims == 1 && !binary) {
            while (getline(f,line)) {
                if (line == "") { continue; };
                iss = std::istringstream(line);
                iss >> x;
                s = _look_up_unit(layer,atoi(x.c_str()));
                for (auto i=0; i<_species_possible_vec[layer].size(); i++) {
                    iss >> sp >> prob;
                    prob_val = atof(prob.c_str());
                    occs_to_write_prob[sp].push_back(std::make_pair(s,prob_val));
                    // Reset
                    sp=""; prob="";
                };
                // Reset
                x="";
            };
        } else if (_no_dims == 2 && binary) {
            while (getline(f,line)) {
                if (line == "") { continue; };
                iss = std::istringstream(line);
                iss >> x >> y >> sp;
                s = _look_up_unit(layer,atoi(x.c_str()),atoi(y.c_str()));
                occs_to_write_binary[sp].push_back(s);
                // Reset
                sp=""; x=""; y="";
            };
        } else if (_no_dims == 2 && !binary) {
            while (getline(f,line)) {
                if (line == "") { continue; };
                iss = std::istringstream(line);
                iss >> x >> y;
                s = _look_up_unit(layer,atoi(x.c_str()),atoi(y.c_str()));
                for (auto i=0; i<_species_possible_vec[layer].size(); i++) {
                    iss >> sp >> prob;
                    prob_val = atof(prob.c_str());
                    occs_to_write_prob[sp].push_back(std::make_pair(s,prob_val));
                    // Reset
                    sp=""; prob="";
                };
                // Reset
                x=""; y="";
            };
        } else if (_no_dims == 3 && binary) {
            while (getline(f,line)) {
                if (line == "") { continue; };
                iss = std::istringstream(line);
                iss >> x >> y >> z >> sp;
                s = _look_up_unit(layer,atoi(x.c_str()),atoi(y.c_str()),atoi(z.c_str()));
                occs_to_write_binary[sp].push_back(s);
                // Reset
                sp=""; x=""; y=""; z="";
            };
        } else if (_no_dims == 3 && !binary) {
            while (getline(f,line)) {
                if (line == "") { continue; };
                iss = std::istringstream(line);
                iss >> x >> y >> z;
                s = _look_up_unit(layer,atoi(x.c_str()),atoi(y.c_str()),atoi(z.c_str()));
                for (auto i=0; i<_species_possible_vec[layer].size(); i++) {
                    iss >> sp >> prob;
                    prob_val = atof(prob.c_str());
                    occs_to_write_prob[sp].push_back(std::make_pair(s,prob_val));
                    // Reset
                    sp=""; prob="";
                };
                // Reset
                x=""; y=""; z="";
            };
        };

        if (binary) {
            // Binary mode
            for (auto const &pr: occs_to_write_binary) {
                // find the species
                auto it = _species_possible_map[layer].find(pr.first);
                if (it == _species_possible_map[layer].end()) {
                    std::cerr << ">>> Lattice::read_from_file <<< Error: could not find species (binary): " << pr.first << std::endl;
                    exit(EXIT_FAILURE);
                };
                for (auto const &site: pr.second) {
                    _mc_chains[chain][i_chain][layer][it->second][site] = 1.0;
                };
            };
        } else {
            // Prob mode
            for (auto const &pr: occs_to_write_prob) {
                // find the species
                auto it = _species_possible_map[layer].find(pr.first);
                if (it == _species_possible_map[layer].end()) {
                    std::cerr << ">>> Lattice::read_from_file <<< Error: could not find species (prob): " << pr.first << std::endl;
                    exit(EXIT_FAILURE);
                };
                for (auto const &site_pr: pr.second) {
                    _mc_chains[chain][i_chain][layer][it->second][site_pr.first] = site_pr.second;
                };
            };
        };
        
		f.close();
	};

    // ***************
    // MARK: - Activate layer steps - PRIVATE
    // ***************
    
    // Reset activations to 0
    void Lattice::_reset_activations(MCType chain, int i_chain, int layer) {
        for (auto sp: _species_possible_vec.at(layer)) {
            _mc_chains_act[chain][i_chain][layer][sp].fill(0.0);
        };
    };
    
    // Calculate activation given layer above or below
    void Lattice::_calculate_activations_from_below(MCType chain, int i_chain, int layer) {
        
        // Get no units
        int no_units = get_no_units_in_layer(layer);
        
        // Activate from below
        for (auto sp: _species_possible_vec.at(layer)) {
            
            // bias term
            _mc_chains_act[chain][i_chain][layer][sp] += get_bias_in_layer(layer, sp) * arma::vec(no_units,arma::fill::ones);
            
            // ixns
            for (auto &given_sp: _species_possible_vec.at(layer-1)) {
                
                // Activate from below
                _mc_chains_act[chain][i_chain][layer][sp] += get_ixn_between_layers(layer-1, given_sp, layer, sp) * ( _adj[layer-1] * _mc_chains[chain][i_chain][layer-1][given_sp] );
            };
        };
    };
    
    void Lattice::_calculate_activations_from_above(MCType chain, int i_chain, int layer) {

        // Get no units
        int no_units = get_no_units_in_layer(layer);

        // Activate from above
        for (auto sp: _species_possible_vec.at(layer)) {
            
            // bias term
            _mc_chains_act[chain][i_chain][layer][sp] += get_bias_in_layer(layer, sp) * arma::vec(no_units,arma::fill::ones);
            
            // ixns
            for (auto given_sp: _species_possible_vec.at(layer+1)) {
                
                // Activate from above
                _mc_chains_act[chain][i_chain][layer][sp] += get_ixn_between_layers(layer+1, given_sp, layer, sp) * ( _adj[layer].t() * _mc_chains[chain][i_chain][layer+1][given_sp] );
            };
        };
    };

    // Calculate activation given layer above or below
    void Lattice::_calculate_activations_from_below_bn(MCType chain, int i_chain, int layer) {
        
        // Get no units
        int no_units = get_no_units_in_layer(layer);
        
        // Activate from below and use batch norm
        if (layer == 1) {
            
            // First hidden layer = eqn (7), part 1, terms 1 and 2
            
            for (auto sp: _species_possible_vec.at(layer)) {
                
                // bias term
                _mc_chains_act[chain][i_chain][layer][sp] += get_bias_in_layer(layer, sp) * arma::vec(no_units,arma::fill::ones);
                
                // ixns
                for (auto &given_sp: _species_possible_vec.at(layer-1)) {
                    
                    // Activate from below
                    _mc_chains_act[chain][i_chain][layer][sp] += get_ixn_between_layers(layer-1, given_sp, layer, sp) * ( _adj[layer-1] * _mc_chains[chain][i_chain][layer-1][given_sp] );
                };
            };
            
        } else {
            
            // Not the first layer = eqn (7), part 2, both terms
            
            for (auto sp: _species_possible_vec.at(layer)) {
                
                // bias term
                _mc_chains_act[chain][i_chain][layer][sp] += get_bias_in_layer(layer, sp) * arma::vec(no_units,arma::fill::ones);
                
                // ixns
                for (auto given_sp: _species_possible_vec.at(layer-1)) {
                    
                    // Activate from below
                    _mc_chains_act[chain][i_chain][layer][sp] += _bn_gamma_bar[chain][layer-1][given_sp] * get_ixn_between_layers(layer-1, given_sp, layer, sp) * ( _adj[layer-1] * _mc_chains[chain][i_chain][layer-1][given_sp] );
                };
            };
        };
    };
    void Lattice::_calculate_activations_from_above_bn(MCType chain, int i_chain, int layer) {
        
        // Get no units
        int no_units = get_no_units_in_layer(layer);
        
        // Activate from above and use batch norm
        
        // Eqn (7), part 1, terms 1 and 3
        // Note: this is equivalent to eqn (4)
        
        for (auto sp: _species_possible_vec.at(layer)) {
            
            // bias term
            _mc_chains_act[chain][i_chain][layer][sp] += get_bias_in_layer(layer, sp) * arma::vec(no_units,arma::fill::ones);
            
            // ixns
            for (auto given_sp: _species_possible_vec.at(layer+1)) {
                
                // Activate from above
                _mc_chains_act[chain][i_chain][layer][sp] += _bn_gamma_bar[chain][layer+1][given_sp] * get_ixn_between_layers(layer+1, given_sp, layer, sp) * ( _adj[layer].t() * _mc_chains[chain][i_chain][layer+1][given_sp] );
            };
        };
    };
    
    // ***************
    // MARK: - Activate layer steps - PUBLIC
    // ***************
    
    // (1) Calculate activations for a specific layer
    // Both directions
    void Lattice::activate_layer_calculate(MCType chain, int i_chain, int layer) {
        _reset_activations(chain,i_chain,layer);
        if (layer != 0) {
            _calculate_activations_from_below(chain,i_chain,layer);
            // _calculate_activations_from_below_bn(chain,i_chain,layer);
        };
        if (layer != _no_layers-1) {
            _calculate_activations_from_above(chain,i_chain,layer);
            // _calculate_activations_from_above_bn(chain,i_chain,layer);
        };
    };
    // Only one direction
    void Lattice::activate_layer_calculate(MCType chain, int i_chain, int layer, int given_layer) {
        _reset_activations(chain,i_chain,layer);
        if (given_layer == layer - 1) {
            _calculate_activations_from_below(chain,i_chain,layer);
            // _calculate_activations_from_below_bn(chain,i_chain,layer);
        } else if (given_layer == layer + 1) {
            _calculate_activations_from_above(chain,i_chain,layer);
            // _calculate_activations_from_above_bn(chain,i_chain,layer);
        } else {
            std::cerr << ">>> Lattice::activate_layer_prepare <<< given layer must be +- layer, but instead layer = " << layer << " and given layer = " << given_layer << std::endl;
            exit(EXIT_FAILURE);
        };
    };
    
    // (2) Convert activations to probs
    void Lattice::activate_layer_convert_to_probs(MCType chain, int i_chain, int layer, bool binary) {
        
        // Convert activations to propensities via exp
        // Also calculate total propensity
        // Starts at 1.0 = exp(0) for empty
        int no_units = get_no_units_in_layer(layer);
        auto prop_tot = arma::vec(no_units,arma::fill::ones);
        for (auto sp: _species_possible_vec.at(layer)) {
            _mc_chains_act[chain][i_chain][layer][sp] = exp(_mc_chains_act[chain][i_chain][layer][sp]);
            prop_tot += _mc_chains_act[chain][i_chain][layer][sp];
        };
        
        // Divide by total to normalize
        // std::cout << "_convert_activations_to_probs: layer: " << layer << " binary: " << binary << std::endl;
        /*
        if (chain == MCType::AWAKE) {
            std::cout << "Lattice::activate_layer_convert_to_probs: awake " << i_chain << " " << layer << " " << binary << std::endl;
        } else {
            std::cout << "Lattice::activate_layer_convert_to_probs: asleep " << i_chain << " " << layer << " " << binary << std::endl;
        };
         */
        for (auto sp: _species_possible_vec.at(layer)) {
            _mc_chains_act[chain][i_chain][layer][sp] /= prop_tot;
            // std::cout << "prob: " << sp->get_name() << " : " << arma::sum(_mc_chains_act[chain][i_chain][layer][sp]) << std::endl;
        };
        
        // Sample if binary
        if (binary) {
            _binarize_all_units_in_layer(chain,i_chain,layer,true);
            /*
            for (auto sp: _species_possible_vec.at(layer)) {
                std::cout << "binary " << sp->get_name() << " : " << arma::sum(_mc_chains_act[chain][i_chain][layer][sp]) << std::endl;
            };
             */
        };
    };
    
    // (3) Commit the new probabilities
    void Lattice::activate_layer_committ(MCType chain, int i_chain, int layer) {
        /*
        std::cout << "activate_layer_committ: layer: " << layer << " current size: " << std::endl;
        for (auto &x: _latt[_i_markov_chain][layer]) {
            std::cout << arma::size(x.second) << std::endl;
        };
        std::cout << "new size:" << std::endl;
        for (auto &x: _latt_act[_i_markov_chain][layer]) {
            std::cout << arma::size(x.second) << std::endl;
        };
         */
        
        for (auto sp: _species_possible_vec.at(layer)) {
            _mc_chains[chain][i_chain][layer][sp] = arma::vec(_mc_chains_act[chain][i_chain][layer][sp]);
        };
    };
    
    // ***************
    // MARK: - Mean field / gibbs sampling
    // ***************
    
    // Activate a specific layer
    void Lattice::activate_single_layer(MCType chain, int i_chain, int layer, bool binary) {
        activate_layer_calculate(chain, i_chain, layer);
        activate_layer_convert_to_probs(chain, i_chain, layer, binary);
        activate_layer_committ(chain, i_chain, layer);
    };
    void Lattice::activate_single_layer(MCType chain, int i_chain, int layer, int given_layer, bool binary) {
        activate_layer_calculate(chain, i_chain, layer, given_layer);
        activate_layer_convert_to_probs(chain, i_chain, layer, binary);
        activate_layer_committ(chain, i_chain, layer);
    };
    
    // Variational inference
    void Lattice::mean_field_hiddens_step() {
        
        // Calculate activations for all chains, all layers
        for (auto i_chain=0; i_chain<_no_markov_chains[MCType::AWAKE]; i_chain++) {
            for (auto layer=1; layer<_no_layers; layer++) {
                activate_layer_calculate(MCType::AWAKE, i_chain, layer);
            };
        };
        
        // Batch normalization
        if (_bn_mode) {
            
            for (auto layer=1; layer<_no_layers; layer++) { // Only for hidden layers
                
                // Calculate means and variances
                _bn_calculate_means_vars_from_activations(MCType::AWAKE, layer);
                
                // Calculate bar params
                _bn_calculate_bar_params_from_means_vars(MCType::AWAKE, layer);
                
                // Apply affine transformations
                _bn_apply_affine_transform_to_all_chains(MCType::AWAKE, layer);
                
            };
        };

        // Convert activations to probabilities
        // Use prob units!
        for (auto i_chain=0; i_chain<_no_markov_chains[MCType::AWAKE]; i_chain++) {
            for (auto layer=1; layer<_no_layers; layer++) {
                activate_layer_convert_to_probs(MCType::AWAKE, i_chain, layer, false);
            };
        };
        
        // Committ new probabilities
        for (auto i_chain=0; i_chain<_no_markov_chains[MCType::AWAKE]; i_chain++) {
            for (auto layer=1; layer<_no_layers; layer++) {
                activate_layer_committ(MCType::AWAKE, i_chain, layer);
            };
        };
    };
    
    // Sample
    void Lattice::gibbs_sampling_step(bool binary_visible, bool binary_hidden) {
        
        // Activate hiddens from visibles
        
        // Calculate activations for all chains, hidden layers
        for (auto i_chain=0; i_chain<_no_markov_chains[MCType::ASLEEP]; i_chain++) {
            for (auto layer=1; layer<_no_layers; layer++) {
                activate_layer_calculate(MCType::ASLEEP, i_chain, layer);
            };
        };

        // Batch normalization
        if (_bn_mode) {
            
            for (auto layer=1; layer<_no_layers; layer++) { // Only for hidden layers
                
                // Calculate means and variances
                _bn_calculate_means_vars_from_activations(MCType::ASLEEP, layer);
                
                // Calculate bar params
                _bn_calculate_bar_params_from_means_vars(MCType::ASLEEP, layer);
                
                // Apply affine transformations
                _bn_apply_affine_transform_to_all_chains(MCType::ASLEEP, layer);
                
            };
        };

        // Convert activations to probabilities
        for (auto i_chain=0; i_chain<_no_markov_chains[MCType::ASLEEP]; i_chain++) {
            for (auto layer=1; layer<_no_layers; layer++) {
                activate_layer_convert_to_probs(MCType::ASLEEP, i_chain, layer, binary_hidden);
            };
        };
        
        // Committ new probabilities
        for (auto i_chain=0; i_chain<_no_markov_chains[MCType::ASLEEP]; i_chain++) {
            for (auto layer=1; layer<_no_layers; layer++) {
                activate_layer_committ(MCType::ASLEEP, i_chain, layer);
            };
        };
        
        // Activate visibles from hiddens
        
        // Calculate activations for all chains, visible layer
        for (auto i_chain=0; i_chain<_no_markov_chains[MCType::ASLEEP]; i_chain++) {
            activate_layer_calculate(MCType::ASLEEP, i_chain, 0);
        };
        
        // Batch normalization is only for hidden layers!
        
        // Convert activations to probabilities
        for (auto i_chain=0; i_chain<_no_markov_chains[MCType::ASLEEP]; i_chain++) {
            activate_layer_convert_to_probs(MCType::ASLEEP, i_chain, 0, binary_visible);
        };
        
        // Committ new probabilities
        for (auto i_chain=0; i_chain<_no_markov_chains[MCType::ASLEEP]; i_chain++) {
            activate_layer_committ(MCType::ASLEEP, i_chain, 0);
        };
    };
    
    // Make a pass activating upwards
    void Lattice::activate_upward_pass(MCType chain, int i_chain, bool binary_hidden) {
        for (auto layer=1; layer<_no_layers; layer++) {
            activate_single_layer(chain, i_chain, layer, layer-1, binary_hidden);
        };
    };
         
    /********************
	Get counts
	********************/
    
    /*
	// 1 particle
	double Lattice::get_count_vis(Sptr &sp) const {
        return sum(_latt.at(_i_markov_chain).at(0).at(sp));
	};

	// 2 particle
	double Lattice::get_count_vis(Sptr &sp1, Sptr &sp2, bool reversibly) const {
		// Only dim=1 for now
		if (_no_dims != 1) {
			std::cerr << ">>> Error: Lattice::get_count <<< only supported for dim 1." << std::endl;
			exit(EXIT_FAILURE);
		};

		UnitVisible *nbr = nullptr;
		double count = 0.0;
		for (auto &s: _latt_v) {
			// Only connect to "plus one" (not minus one) => no duplicates!

			if (s->x()+1 <= _box_length) {
				// Get neighbor
				nbr = _look_up_unit_v(s->x()+1);

				// Count
				count += s->get_occ(sp1) * nbr->get_occ(sp2); 
				if (reversibly && sp1 != sp2) {
					count += nbr->get_occ(sp1) * s->get_occ(sp2); 
				};
			};
		};

		return count;
	};

	// 3 particle
	double Lattice::get_count_vis(Sptr &sp1, Sptr &sp2, Sptr &sp3, bool reversibly) const {
		// Only dim=1 for now
		if (_no_dims != 1) {
			std::cerr << ">>> Error: Lattice::get_count <<< only supported for dim 1." << std::endl;
			exit(EXIT_FAILURE);
		};

		UnitVisible *nbr1 = nullptr, *nbr2 = nullptr;
		double count = 0.0;
		for (auto &s: _latt_v) {
			// Only connect to "plus one" (not minus one) => no duplicates!

			if (s->x()+2 <= _box_length) {
				// Get neighbors
				nbr1 = _look_up_unit_v(s->x()+1);
				nbr2 = _look_up_unit_v(s->x()+2);

				// Count
				count += s->get_occ(sp1) * nbr1->get_occ(sp2) * nbr2->get_occ(sp3); 
				if (reversibly && sp1 != sp3) {
					count += nbr2->get_occ(sp1) * nbr1->get_occ(sp2) * s->get_occ(sp3); 
				};
			};
		};

		return count;
	};

	// 4 particle
	double Lattice::get_count_vis(Sptr &sp1, Sptr &sp2, Sptr &sp3, Sptr &sp4, bool reversibly) const {
		// Only dim=1 for now
		if (_no_dims != 1) {
			std::cerr << ">>> Error: Lattice::get_count <<< only supported for dim 1." << std::endl;
			exit(EXIT_FAILURE);
		};

		UnitVisible *nbr1 = nullptr, *nbr2 = nullptr, *nbr3 = nullptr;
		double count = 0.0;
		for (auto &s: _latt_v) {
			// Only connect to "plus one" (not minus one) => no duplicates!

			if (s->x()+3 <= _box_length) {
				// Get neighbors
				nbr1 = _look_up_unit_v(s->x()+1);
				nbr2 = _look_up_unit_v(s->x()+2);
				nbr3 = _look_up_unit_v(s->x()+3);

				// Count
				count += s->get_occ(sp1) * nbr1->get_occ(sp2) * nbr2->get_occ(sp3) * nbr3->get_occ(sp4); 
				if (reversibly && !(sp1 == sp4 && sp2 == sp3)) {
					count += nbr3->get_occ(sp1) * nbr2->get_occ(sp2) * nbr1->get_occ(sp3) * s->get_occ(sp4); 
				};
			};
		};

		return count;
	};
     */
    
    // ***************
    // MARK: - Reap moments
    // ***************
    
    void Lattice::reap_moments(MCType chain) {
        // Reset all ixns
        for (auto ixn: _all_ixns) {
            for (auto i_chain=0; i_chain<_no_markov_chains.at(chain); i_chain++) {
                ixn->get_moment()->set_moment_sample(chain, i_chain, 0.0);
            };
        };
        
        // Reap biases
        double val;
        // Go over all chains
        for (auto i_chain=0; i_chain<_no_markov_chains.at(chain); i_chain++) {
            // All biases
            for (auto &bias_layer: _bias_dict) {
                // Go through possible species in this layer
                for (auto &sp_pr: bias_layer.second) {
                    // Get the moment
                    val = arma::sum(_mc_chains.at(chain).at(i_chain).at(bias_layer.first).at(sp_pr.first));
            
                    // Go through all ixns associated with this species
                    for (auto ixn: sp_pr.second) {
                        // Set the moment
                        if (ixn->get_moment()) {
                            if (chain == MCType::ASLEEP || !ixn->get_moment()->get_is_awake_moment_fixed()) {
                                ixn->get_moment()->increment_moment_sample(chain, i_chain, val);
                            };
                        };
                    };
                };
            };
        };
        
        // Reap ixns
        // Go over all chains
        for (auto i_chain=0; i_chain<_no_markov_chains.at(chain); i_chain++) {
            // All ixns
            for (auto &o2_ixn_layer_1: _o2_ixn_dict) {
                // Go through possible species in this layer
                for (auto &sp_pr_1: o2_ixn_layer_1.second) {
                    for (auto &o2_ixn_layer_2: sp_pr_1.second) {
                        // Be careful not to double count
                        if (o2_ixn_layer_2.first <= o2_ixn_layer_1.first) {
                            // skip
                            continue;
                        };
                        
                        // Now ( layer 1 = o2_ixn_layer_1.first ) < ( layer 2 = o2_ixn_layer_2.first ) guaranteed
                        
                        for (auto &sp_pr_2: o2_ixn_layer_2.second) {
                            // Get the moment
                            val = dot(_mc_chains.at(chain).at(i_chain).at(o2_ixn_layer_1.first).at(sp_pr_1.first), _adj.at(o2_ixn_layer_1.first) * _mc_chains.at(chain).at(i_chain).at(o2_ixn_layer_2.first).at(sp_pr_2.first));
                            
                            // Go through all ixns associated with this species
                            for (auto ixn: sp_pr_2.second) {
                                // Set the moment
                                if (ixn->get_moment()) {
                                    if (chain == MCType::ASLEEP || !ixn->get_moment()->get_is_awake_moment_fixed()) {
                                        ixn->get_moment()->increment_moment_sample(chain, i_chain, val);
                                    };
                                };
                            };
                        };
                    };
                };
            };
        };
        
        // Gamma, beta
        if (_bn_mode) {
            
            // Activate all chains again first (but not the visible layer!)
            for (auto layer=1; layer<_no_layers; layer++) {
                for (auto i_chain=0; i_chain<_no_markov_chains.at(chain); i_chain++) {
                    activate_layer_calculate(chain, i_chain, layer);
                    /*
                    if (chain == MCType::AWAKE) {
                        for (auto sp: _species_possible_vec.at(layer)) {
                            std::cout << "activated init: " << layer << " " << i_chain << " " << arma::sum(_mc_chains_act.at(chain).at(i_chain).at(layer).at(sp)) << std::endl;
                        };
                    };
                     */
                };
            };
            
            // Calculate mean, var
            for (auto layer=1; layer<_no_layers; layer++) {
                _bn_calculate_means_vars_from_activations(chain, layer);
                
                // Transform all chains
                for (auto i_chain=0; i_chain<_no_markov_chains.at(chain); i_chain++) {
                    for (auto sp: _species_possible_vec[layer]) {
                        
                        /*
                        if (chain == MCType::AWAKE) {
                            std::cout << "initial: " << arma::sum(_mc_chains_act.at(chain).at(i_chain).at(layer).at(sp)) << std::endl;
                            std::cout << "means: " << arma::sum(_bn_means[chain][layer][sp]) << std::endl;
                        };
                         */
                        
                        _mc_chains_act[chain][i_chain][layer][sp] -= _bn_means[chain][layer][sp];
                        
                        /*
                        if (chain == MCType::AWAKE) {
                            std::cout << "subtracted: " << arma::sum(_mc_chains_act.at(chain).at(i_chain).at(layer).at(sp)) << std::endl;
                        };
                         */
                        
                        _mc_chains_act[chain][i_chain][layer][sp] /= sqrt(_bn_vars[chain][layer][sp] + _bn_eps);
                        
                        /*
                        if (chain == MCType::AWAKE) {
                            std::cout << "divided: " << arma::sum(_mc_chains_act.at(chain).at(i_chain).at(layer).at(sp)) << std::endl;
                        };
                         */
                    };
                };
            };
            
            // Assign
            double val_gamma,val_beta;
            for (auto layer=1; layer<_no_layers; layer++) {
                for (auto sp: _species_possible_vec.at(layer)) {
                    for (auto i_chain=0; i_chain<_no_markov_chains.at(chain); i_chain++) {
                        // Calculate
                        /*
                        if (chain == MCType::AWAKE) {
                            std::cout << "activated: " << arma::sum(_mc_chains_act.at(chain).at(i_chain).at(layer).at(sp)) << std::endl;
                        };
                         */
                        val_gamma = arma::dot(_mc_chains.at(chain).at(i_chain).at(layer).at(sp), _mc_chains_act.at(chain).at(i_chain).at(layer).at(sp));
                        val_beta = arma::sum(_mc_chains.at(chain).at(i_chain).at(layer).at(sp));
                        // Assign
                        _bn_gamma[layer][sp]->get_moment()->set_moment_sample(chain, i_chain, val_gamma);
                        _bn_beta[layer][sp]->get_moment()->set_moment_sample(chain, i_chain, val_beta);
                    };
                };
            };
        };
    };
    
    // *******************
    // MARK: - Batch normalization
    // *******************
    
    // Calculate means from activations
    void Lattice::_bn_calculate_means_vars_from_activations(MCType chain, int layer) {
        // Reset
        for (auto sp: _species_possible_vec[layer]) {
            _bn_means[chain][layer][sp].fill(0.0);
            _bn_vars[chain][layer][sp].fill(0.0);
        };
        
        // Mean
        for (auto sp: _species_possible_vec[layer]) {
            for (auto i_chain=0; i_chain<_no_markov_chains[chain]; i_chain++) {
                _bn_means[chain][layer][sp] += _mc_chains_act[chain][i_chain][layer][sp];
            };
            _bn_means[chain][layer][sp] /= _no_markov_chains[chain];
        };
        
        // Variance
        for (auto sp: _species_possible_vec[layer]) {
            for (auto i_chain=0; i_chain<_no_markov_chains[chain]; i_chain++) {
                _bn_vars[chain][layer][sp] += pow(_mc_chains_act[chain][i_chain][layer][sp] - _bn_means[chain][layer][sp],2);
            };
            _bn_vars[chain][layer][sp] /= _no_markov_chains[chain];
        };
    };
    
    // Calculate bar parameters from means, vars
    void Lattice::_bn_calculate_bar_params_from_means_vars(MCType chain, int layer) {
        for (auto sp: _species_possible_vec[layer]) {
            _bn_beta_bar[chain][layer][sp] = _bn_beta[layer][sp]->get_val() - ( _bn_gamma[layer][sp]->get_val() * _bn_means[chain][layer][sp] / sqrt(_bn_vars[chain][layer][sp] + _bn_eps) );
            _bn_gamma_bar[chain][layer][sp] = _bn_gamma[layer][sp]->get_val() / sqrt(_bn_vars[chain][layer][sp] + _bn_eps);
        };
    };

    // Apply BN transform
    void Lattice::_bn_apply_affine_transform_to_all_chains(MCType chain, int layer) {
        for (auto i_chain=0; i_chain<_no_markov_chains[chain]; i_chain++) {
            for (auto sp: _species_possible_vec[layer]) {
                _mc_chains_act[chain][i_chain][layer][sp] %= _bn_gamma_bar[chain][layer][sp];
                _mc_chains_act[chain][i_chain][layer][sp] += _bn_beta_bar[chain][layer][sp];
            };
        };
    };
};
