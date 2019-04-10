#include "../include/dblz_bits/lattice.hpp"

// Other headers
#include "../include/dblz_bits/general.hpp"
#include "../include/dblz_bits/species.hpp"
#include "../include/dblz_bits/moment.hpp"
#include "../include/dblz_bits/ixn_param.hpp"
#include "../include/dblz_bits/fname.hpp"

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
	Lattice
	****************************************/

	/********************
	Constructor
	********************/

    Lattice::Lattice(int no_dims, int box_length, std::vector<Sptr> species_visible, LatticeMode mode) {
        
        if (no_dims != 1 && no_dims != 2 && no_dims != 3) {
            std::cerr << "ERROR: only dimensions 1,2,3 are supported for Lattice." << std::endl;
            exit(EXIT_FAILURE);
        };
        
        _no_dims = no_dims;
        _box_length = box_length;
        _no_layers = 0;
        _mode = mode;
        
        // Set no chains
        set_no_markov_chains(MCType::AWAKE, 1);
        set_no_markov_chains(MCType::ASLEEP, 1);
        
        // Visible layer
        add_layer(0, _box_length, species_visible);
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
        _bias_mults = other._bias_mults;
        
        _mode = other._mode;
        
        _c_means = other._c_means;
        _c_batch_means = other._c_batch_means;
        
        _cpt_means = other._cpt_means;
        _cpt_batch_means = other._cpt_batch_means;
        
        _pst_prop = other._pst_prop;
        _pst_r = other._pst_r;
        _pst_sign_of_r = other._pst_sign_of_r;
        _pst_sign_of_r_new = other._pst_sign_of_r_new;
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
        _bias_mults = other._bias_mults;

        _mode = other._mode;
        
        _c_means = other._c_means;
        _c_batch_means = other._c_batch_means;
        
        _cpt_means = other._cpt_means;
        _cpt_batch_means = other._cpt_batch_means;
    
        _pst_prop = other._pst_prop;
        _pst_r = other._pst_r;
        _pst_sign_of_r = other._pst_sign_of_r;
        _pst_sign_of_r_new = other._pst_sign_of_r_new;
        
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
        other._bias_mults.clear();
        
        other._c_means.clear();
        other._c_batch_means.clear();
        
        other._cpt_means.clear();
        other._cpt_batch_means.clear();
        
        other._pst_prop.clear();
        other._pst_r.clear();
        other._pst_sign_of_r.clear();
        other._pst_sign_of_r_new.clear();
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
    
    LatticeMode Lattice::get_lattice_mode() const {
        return _mode;
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
        int ctr=0;
        if (_no_dims == 1) {
            for (auto x=1; x<=box_length; x++) {
                _lookup_1[layer][x] = ctr;
                _rlookup[layer][ctr] = std::vector<int>({x});
                ctr++;
            };
        } else if (_no_dims == 2) {
            for (auto x=1; x<=box_length; x++) {
                for (auto y=1; y<=box_length; y++) {
                    _lookup_2[layer][x][y] = ctr;
                    _rlookup[layer][ctr] = std::vector<int>({x,y});
                    ctr++;
                };
            };
        } else if (_no_dims == 3) {
            for (auto x=1; x<=box_length; x++) {
                for (auto y=1; y<=box_length; y++) {
                    for (auto z=1; z<=box_length; z++) {
                        // std::cout << "Setting idx in layer: " << layer << " pos: " << x << " " << y << " " << z << " to idx: " << ctr << std::endl;
                        _lookup_3[layer][x][y][z] = ctr;
                        _rlookup[layer][ctr] = std::vector<int>({x,y,z});
                        ctr++;
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
            _adj[layer-1][layer] = arma::sp_mat(no_units,size_below);
            _adj[layer][layer-1] = _adj[layer-1][layer].t();
            // std::cout << "Made adjacency matrix: " << layer-1 << " " << no_units << " " << size_below << std::endl;
        };

        // Centered?
        if (_mode == LatticeMode::NORMAL_W_CENTERED_GRADIENT_VEC) {
            int no_units = get_no_units_in_layer(layer);
            for (auto sp: _species_possible_vec.at(layer)) {
                _c_means[layer][sp] = arma::vec(no_units,arma::fill::ones);
                _c_batch_means[layer][sp] = arma::vec(no_units,arma::fill::ones);
                
                // Init to 0.5
                _c_means[layer][sp] *= 0.5;
                _c_batch_means[layer][sp] *= 0.5;
            };
        } else if (_mode == LatticeMode::NORMAL_W_CENTERED_GRADIENT_PT || _mode == LatticeMode::CENTERED_PT) {
            for (auto sp: _species_possible_vec.at(layer)) {
                _cpt_means[layer][sp] = 0.5;
                _cpt_batch_means[layer][sp] = 0.5;
            };
        };
        
        // Persistent data structures
        _pst_prop[layer] = arma::vec(no_units);
        _pst_r[layer] = arma::vec(no_units);
        _pst_sign_of_r[layer] = arma::vec(no_units);
        _pst_sign_of_r_new[layer] = arma::vec(no_units);
    };
    
	/********************
	Helpers to setup all sites - Biases
	********************/

    // Biases
    void Lattice::set_bias_of_layer(int layer, Sptr sp, Iptr bias) {
        _bias_dict[layer][sp] = bias;
        
        _add_to_all_ixns_vec(bias);
        
        if (_mode == LatticeMode::NORMAL_W_CENTERED_GRADIENT_VEC) {
            int no_units = get_no_units_in_layer(layer);
            bias->get_moment()->set_bias_vec_dims(no_units);
        };
    };

    // Ixns
    void Lattice::set_ixn_between_layers(int layer1, Sptr sp1, int layer2, Sptr sp2, Iptr ixn) {
        if (layer2 != layer1+1 && layer2 != layer1-1) {
            std::cerr << ">>> Lattice::add_ixn_between_layers <<< layer2 != layer1 +- 1; instead layer_2 = " << layer2 << " and layer1 = " << layer1 << std::endl;
            exit(EXIT_FAILURE);
        };
        
        // Add both ways
        _o2_ixn_dict[layer1][sp1][layer2][sp2] = ixn;
        _o2_ixn_dict[layer2][sp2][layer1][sp1] = ixn;

        // Add to all
        _add_to_all_ixns_vec(ixn);
        
        // Set moment matrix dims if we are in centered mode
        if (_mode == LatticeMode::NORMAL_W_CENTERED_GRADIENT_VEC) {
            int no_below = get_no_units_in_layer(std::min(layer1,layer2));
            int no_above = get_no_units_in_layer(std::max(layer1,layer2));
            ixn->get_moment()->set_weight_matrix_dims(no_above,no_below);
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
    void Lattice::set_multiplier_for_bias_in_layer(int layer, double multiplier) {
        _bias_mults[layer] = multiplier;
    };
    
    // Get ixns
    double Lattice::get_bias_in_layer(int layer, Sptr sp) const {
        
        // Multiplier
        double mult = 1.0;
        auto itm = _bias_mults.find(layer);
        if (itm != _bias_mults.end()) {
            mult = itm->second;
        };
        
        auto it = _bias_dict.find(layer);
        double val = 0.0;
        if (it != _bias_dict.end()) {
            auto it2 = it->second.find(sp);
            if (it2 != it->second.end()) {
                val = it2->second->get_val();
            };
        };
        return mult * val;
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
                        val = it4->second->get_val();
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
        _adj[layer1][layer2](idx2,idx1) = 1.0;
        _adj[layer2][layer1](idx1,idx2) = 1.0;
    };
    void Lattice::add_conn(int layer1, int x1, int y1, int layer2, int x2, int y2) {
        int idx1 = _look_up_unit(layer1, x1, y1);
        int idx2 = _look_up_unit(layer2, x2, y2);
        _adj[layer1][layer2](idx2,idx1) = 1.0;
        _adj[layer2][layer1](idx1,idx2) = 1.0;
    };
    void Lattice::add_conn(int layer1, int x1, int y1, int z1, int layer2, int x2, int y2, int z2) {
        int idx1 = _look_up_unit(layer1, x1, y1, z1);
        int idx2 = _look_up_unit(layer2, x2, y2, z2);
        // std::cout << "Connecting: " << layer1 << " " << x1 << " " << y1 << " " << z1 << " : " << layer2 << " " << x2 << " " << y2 << " " << z2 << " : " << idx1 << " " << idx2 << std::endl;
        _adj[layer1][layer2](idx2,idx1) = 1.0;
        _adj[layer2][layer1](idx1,idx2) = 1.0;
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
        _pst_r.at(layer).fill(arma::fill::randu);
        
        // Helpers
        
        // sign of r vector:
        // 1 initially
        // 1 -> 1 => not this species
        // 1 -> -1 => this species
        // -1 -> -1 => not this species
        _pst_sign_of_r.at(layer).fill(arma::fill::ones);
        _pst_sign_of_r_new.at(layer).fill(arma::fill::ones);

        // _latt or _latt_act?
        if (act) {
            
            // Evaluate
            // Initially: we (flip) are at zero, rand value is somewhere above
            // Next: we (flip) are at some propensity, rand value may be somewhere above or below
            // Finally: we (flip) are at max propensity = 1.0, rand value is below
            for (auto sp: _species_possible_vec.at(layer)) {
                // Subtract prob from r
                _pst_r[layer] -= _mc_chains_act[chain][i_chain][layer][sp];
                
                // New flip vector
                _pst_sign_of_r_new[layer] = sign(_pst_r.at(layer));
                
                // Difference
                // sign_of_r -> sign_of_r_new
                // 1 -> 1 => not this species [0]
                // 1 -> -1 => this species [1]
                // -1 -> -1 => not this species [0]
                // 0.5 * (sign_of_r - sign_of_r_new) encodes the [brackets]
                
                _mc_chains_act[chain][i_chain][layer][sp] = 0.5 * (_pst_sign_of_r.at(layer) - _pst_sign_of_r_new.at(layer));
                
                // Old flip vector = new flip vector
                _pst_sign_of_r[layer] = _pst_sign_of_r_new.at(layer);
            };
            // Finally: if still above the final propensity, it is the empty species! This means all are zero, so no further adjustment needed!
            
        } else {

            for (auto sp: _species_possible_vec.at(layer)) {
                _pst_r[layer] -= _mc_chains[chain][i_chain][layer][sp];
                _pst_sign_of_r_new[layer] = sign(_pst_r.at(layer));
                _mc_chains[chain][i_chain][layer][sp] = 0.5 * (_pst_sign_of_r.at(layer) - _pst_sign_of_r_new.at(layer));
                _pst_sign_of_r[layer] = _pst_sign_of_r_new.at(layer);
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
    // MARK: - Activate layer steps
    // ***************
    
    // Calculate activation given layer above or below
    void Lattice::activate_layer_calculate_from_below(MCType chain, int layer) {
        
        if (_mode == LatticeMode::NORMAL || _mode == LatticeMode::NORMAL_W_CENTERED_GRADIENT_PT || _mode == LatticeMode::NORMAL_W_CENTERED_GRADIENT_VEC) {
            
            for (auto i_chain=0; i_chain<_no_markov_chains.at(chain); i_chain++) {
                for (auto sp: _species_possible_vec.at(layer)) {
                    
                    // bias term
                    _mc_chains_act[chain][i_chain][layer][sp].fill(get_bias_in_layer(layer, sp));
                    
                    // ixns
                    for (auto &given_sp: _species_possible_vec.at(layer-1)) {
                        
                        // Activate from below
                        _mc_chains_act[chain][i_chain][layer][sp] += get_ixn_between_layers(layer-1, given_sp, layer, sp) * ( _adj.at(layer-1).at(layer) * _mc_chains.at(chain).at(i_chain).at(layer-1).at(given_sp) );
                    };
                };
            };

        } else if (_mode == LatticeMode::CENTERED_PT) {
            
            for (auto i_chain=0; i_chain<_no_markov_chains.at(chain); i_chain++) {
                
                for (auto sp: _species_possible_vec.at(layer)) {
                    
                    // bias term
                    _mc_chains_act[chain][i_chain][layer][sp].fill(get_bias_in_layer(layer, sp));
                    
                    // ixns
                    for (auto &given_sp: _species_possible_vec.at(layer-1)) {
                        
                        // Activate from below
                        _mc_chains_act[chain][i_chain][layer][sp] += get_ixn_between_layers(layer-1, given_sp, layer, sp) * ( _adj.at(layer-1).at(layer) * ( _mc_chains.at(chain).at(i_chain).at(layer-1).at(given_sp) - _cpt_means.at(layer-1).at(given_sp) ) );
                    };
                };
            };
        };
    };

    void Lattice::activate_layer_calculate_from_above(MCType chain, int layer) {
        
        if (_mode == LatticeMode::NORMAL || _mode == LatticeMode::NORMAL_W_CENTERED_GRADIENT_PT || _mode == LatticeMode::NORMAL_W_CENTERED_GRADIENT_VEC) {

            for (auto i_chain=0; i_chain<_no_markov_chains.at(chain); i_chain++) {
                for (auto sp: _species_possible_vec.at(layer)) {
                    
                    // bias term
                    _mc_chains_act[chain][i_chain][layer][sp].fill(get_bias_in_layer(layer, sp));
                    
                    // ixns
                    for (auto given_sp: _species_possible_vec.at(layer+1)) {
                        
                        // Activate from above
                        _mc_chains_act[chain][i_chain][layer][sp] += get_ixn_between_layers(layer+1, given_sp, layer, sp) * ( _adj.at(layer+1).at(layer) * _mc_chains.at(chain).at(i_chain).at(layer+1).at(given_sp) );
                    };
                };
            };

        } else if (_mode == LatticeMode::CENTERED_PT) {
            
            for (auto i_chain=0; i_chain<_no_markov_chains.at(chain); i_chain++) {
                
                for (auto sp: _species_possible_vec.at(layer)) {
                    
                    // bias term
                    _mc_chains_act[chain][i_chain][layer][sp].fill(get_bias_in_layer(layer, sp));
                    
                    // ixns
                    for (auto given_sp: _species_possible_vec.at(layer+1)) {
                        
                        // Activate from above
                        _mc_chains_act[chain][i_chain][layer][sp] += get_ixn_between_layers(layer+1, given_sp, layer, sp) * ( _adj.at(layer+1).at(layer) * ( _mc_chains.at(chain).at(i_chain).at(layer+1).at(given_sp) - _cpt_means.at(layer+1).at(given_sp) ) );
                    };
                };
            };
        };
    };

    void Lattice::activate_layer_calculate_from_both(MCType chain, int layer) {
    
        if (_mode == LatticeMode::NORMAL || _mode == LatticeMode::NORMAL_W_CENTERED_GRADIENT_PT || _mode == LatticeMode::NORMAL_W_CENTERED_GRADIENT_VEC) {
            
            for (auto i_chain=0; i_chain<_no_markov_chains.at(chain); i_chain++) {
                for (auto sp: _species_possible_vec.at(layer)) {
                    
                    // bias term
                    _mc_chains_act[chain][i_chain][layer][sp].fill(get_bias_in_layer(layer, sp));
                    
                    // Activate from above
                    for (auto given_sp: _species_possible_vec.at(layer+1)) {
                        _mc_chains_act[chain][i_chain][layer][sp] += get_ixn_between_layers(layer+1, given_sp, layer, sp) * ( _adj.at(layer+1).at(layer) * _mc_chains.at(chain).at(i_chain).at(layer+1).at(given_sp) );
                    };
                    
                    // Activate from below
                    for (auto given_sp: _species_possible_vec.at(layer-1)) {
                        _mc_chains_act[chain][i_chain][layer][sp] += get_ixn_between_layers(layer-1, given_sp, layer, sp) * ( _adj.at(layer-1).at(layer) * _mc_chains.at(chain).at(i_chain).at(layer-1).at(given_sp) );
                    };
                };
            };

        } else if (_mode == LatticeMode::CENTERED_PT) {
            
            for (auto i_chain=0; i_chain<_no_markov_chains.at(chain); i_chain++) {
                for (auto sp: _species_possible_vec.at(layer)) {
                    
                    // bias term
                    _mc_chains_act[chain][i_chain][layer][sp].fill(get_bias_in_layer(layer, sp));
                    
                    // Activate from above
                    for (auto given_sp: _species_possible_vec.at(layer+1)) {
                        _mc_chains_act[chain][i_chain][layer][sp] += get_ixn_between_layers(layer+1, given_sp, layer, sp) * ( _adj.at(layer+1).at(layer) * ( _mc_chains.at(chain).at(i_chain).at(layer+1).at(given_sp) - _cpt_means.at(layer+1).at(given_sp) ) );
                    };
                    
                    // Activate from below
                    for (auto given_sp: _species_possible_vec.at(layer-1)) {
                        _mc_chains_act[chain][i_chain][layer][sp] += get_ixn_between_layers(layer-1, given_sp, layer, sp) * ( _adj.at(layer-1).at(layer) * ( _mc_chains.at(chain).at(i_chain).at(layer-1).at(given_sp) - _cpt_means.at(layer-1).at(given_sp) ) );
                    };
                };
            };
        };
    };

    // (2) Convert activations to probs
    void Lattice::activate_layer_convert_to_probs(MCType chain, int layer, bool binary) {
        
        // Starts at 1.0 = exp(0) for empty
        _pst_prop.at(layer).fill(arma::fill::ones);

        // All chains
        for (auto i_chain=0; i_chain<_no_markov_chains.at(chain); i_chain++) {
        
            // Convert activations to propensities via exp
            // Also calculate total propensity
            _pst_prop.at(layer).fill(arma::fill::ones);
            for (auto sp: _species_possible_vec.at(layer)) {
                _mc_chains_act[chain][i_chain][layer][sp].transform( [](double val) { return (exp(val)); } );
                _pst_prop[layer] += _mc_chains_act[chain][i_chain][layer][sp];
            };
            
            // Divide by total to normalize
            for (auto sp: _species_possible_vec.at(layer)) {
                _mc_chains_act[chain][i_chain][layer][sp] /= _pst_prop.at(layer);
            };
        };
            
        // Sample if binary
        if (binary) {
            for (auto i_chain=0; i_chain<_no_markov_chains.at(chain); i_chain++) {
                _binarize_all_units_in_layer(chain,i_chain,layer,true);
            };
        };
    };
    
    // (3) Commit the new probabilities
    void Lattice::activate_layer_committ(MCType chain, int layer) {
        // All chains
        for (auto i_chain=0; i_chain<_no_markov_chains.at(chain); i_chain++) {
            for (auto sp: _species_possible_vec.at(layer)) {
                _mc_chains[chain][i_chain][layer][sp] = _mc_chains_act[chain][i_chain][layer][sp];
            };
        };
    };
    
    // ***************
    // MARK: - Mean field / gibbs sampling
    // ***************
    
    // Variational inference
    void Lattice::mean_field_hiddens_step() {
        
        // Activate all layers
        for (auto layer=1; layer<_no_layers-1; layer++) {
            activate_layer_calculate_from_both(MCType::AWAKE, layer);
            activate_layer_convert_to_probs(MCType::AWAKE, layer, false);
            activate_layer_committ(MCType::AWAKE, layer);
        };
        // Last layer
        activate_layer_calculate_from_below(MCType::AWAKE, _no_layers-1);
        activate_layer_convert_to_probs(MCType::AWAKE, _no_layers-1, false);
        activate_layer_committ(MCType::AWAKE, _no_layers-1);
    };

    // Sample
    void Lattice::gibbs_sampling_step_awake(bool binary_hidden) {
        
        // Activate in two blocks: odds and evens!
        
        // First the odd layers
        for (auto layer=1; layer<_no_layers; layer += 2) {
            if (layer != _no_layers-1) {
                activate_layer_calculate_from_both(MCType::AWAKE, layer);
            } else {
                activate_layer_calculate_from_below(MCType::AWAKE, layer);
            };
            activate_layer_convert_to_probs(MCType::AWAKE, layer, binary_hidden);
            activate_layer_committ(MCType::AWAKE, layer);
        };
        
        // Next the even layers
        // Other layers
        for (auto layer=2; layer<_no_layers; layer += 2) {
            if (layer != _no_layers-1) {
                activate_layer_calculate_from_both(MCType::AWAKE, layer);
            } else {
                activate_layer_calculate_from_below(MCType::AWAKE, layer);
            };
            activate_layer_convert_to_probs(MCType::AWAKE, layer, binary_hidden);
            activate_layer_committ(MCType::AWAKE, layer);
        };
    };
    
    void Lattice::gibbs_sampling_step_parallel_awake(bool binary_hidden) {
        
        for (auto layer=1; layer<_no_layers-1; layer++) {
            activate_layer_calculate_from_both(MCType::AWAKE, layer);
        };
        activate_layer_calculate_from_below(MCType::AWAKE, _no_layers-1);
        
        // Convert activations to probabilities and committ
        for (auto layer=1; layer<_no_layers; layer++) {
            activate_layer_convert_to_probs(MCType::AWAKE, layer, binary_hidden);
            activate_layer_committ(MCType::AWAKE, layer);
        };
    };

    // Sample
    void Lattice::gibbs_sampling_step(bool binary_visible, bool binary_hidden) {
        
        // Activate in two blocks: odds and evens!
        
        // First the odd layers
        for (auto layer=1; layer<_no_layers; layer += 2) {
            if (layer != _no_layers-1) {
                activate_layer_calculate_from_both(MCType::ASLEEP, layer);
            } else {
                activate_layer_calculate_from_below(MCType::ASLEEP, layer);
            };
            activate_layer_convert_to_probs(MCType::ASLEEP, layer, binary_hidden);
            activate_layer_committ(MCType::ASLEEP, layer);
        };
        
        // Next the even layers
        // Zeroth layer
        activate_layer_calculate_from_above(MCType::ASLEEP, 0);
        activate_layer_convert_to_probs(MCType::ASLEEP, 0, binary_visible);
        activate_layer_committ(MCType::ASLEEP, 0);
        // Other layers
        for (auto layer=2; layer<_no_layers; layer += 2) {
            if (layer != _no_layers-1) {
                activate_layer_calculate_from_both(MCType::ASLEEP, layer);
            } else {
                activate_layer_calculate_from_below(MCType::ASLEEP, layer);
            };
            activate_layer_convert_to_probs(MCType::ASLEEP, layer, binary_hidden);
            activate_layer_committ(MCType::ASLEEP, layer);
        };
    };

    void Lattice::gibbs_sampling_step_parallel(bool binary_visible, bool binary_hidden) {
        
        // Activate in parallel
        
        activate_layer_calculate_from_above(MCType::ASLEEP, 0);
        for (auto layer=1; layer<_no_layers-1; layer++) {
            activate_layer_calculate_from_both(MCType::ASLEEP, layer);
        };
        activate_layer_calculate_from_below(MCType::ASLEEP, _no_layers-1);
        
        // Convert activations to probabilities and committ
        activate_layer_convert_to_probs(MCType::ASLEEP, 0, binary_visible);
        activate_layer_committ(MCType::ASLEEP, 0);
        for (auto layer=1; layer<_no_layers; layer++) {
            activate_layer_convert_to_probs(MCType::ASLEEP, layer, binary_hidden);
            activate_layer_committ(MCType::ASLEEP, layer);
        };
    };

    // Make a pass activating upwards
    void Lattice::activate_upward_pass(MCType chain, bool binary_hidden) {
        
        // All layers
        for (auto layer=1; layer<_no_layers; layer++) {
            activate_layer_calculate_from_below(chain, layer);
            activate_layer_convert_to_probs(chain, layer, binary_hidden);
            activate_layer_committ(chain, layer);
        };
    };

    void Lattice::activate_upward_pass_with_2x_weights_1x_bias(MCType chain, bool binary_hidden) {
        // Copy over the mults for safekeeping
        std::map<int, std::map<int, double>> o2_mults = _o2_mults;
        std::map<int, double> bias_mults = _bias_mults;
        
        // 2x
        for (auto layer=0; layer<_no_layers-2; layer++) {
            _o2_mults[layer][layer+1] = 2.0;
        };
        // 1x last layer
        _o2_mults[_no_layers-2][_no_layers-1] = 1.0;

        // 1x bias
        for (auto layer=0; layer<_no_layers; layer++) {
            _bias_mults[layer] = 1.0;
        };
        
        // Activate upward pass
        activate_upward_pass(chain, binary_hidden);
        
        // Copy back
        _o2_mults = o2_mults;
        _bias_mults = bias_mults;
    };
    
    // ***************
    // MARK: - Reap moments
    // ***************
    
    void Lattice::reap_moments_and_slide_centers_normal() {
        
        // Reap ixns
        int layer1, layer2;
        Sptr sp1, sp2;
        arma::sp_mat::iterator mit, mit_end;
        std::shared_ptr<Moment> moment;
        for (auto &o2_ixn_layer_1: _o2_ixn_dict) {
            layer1 = o2_ixn_layer_1.first;
            for (auto &sp_pr_1: o2_ixn_layer_1.second) {
                sp1 = sp_pr_1.first;
                for (auto &o2_ixn_layer_2: sp_pr_1.second) {
                    layer2 = o2_ixn_layer_2.first;
                    // Be careful not to double count
                    if (layer1 >= layer2) {
                        // skip
                        continue;
                    };
                    // Now layer1 < layer2 guaranteed
                    
                    for (auto &sp_pr_2: o2_ixn_layer_2.second) {
                        
                        sp2 = sp_pr_2.first;
                        moment = sp_pr_2.second->get_moment();
                        
                        // Reset moments
                        if (!moment->get_is_awake_moment_fixed()) {
                            moment->reset_moment(MCType::AWAKE);
                        };
                        moment->reset_moment(MCType::ASLEEP);

                        // Iterate over adj matrix
                        mit = _adj.at(layer1).at(layer2).begin();
                        mit_end = _adj.at(layer1).at(layer2).end();
                        for(; mit != mit_end; ++mit) {
                            
                            // Iterate over awake chains
                            if (!moment->get_is_awake_moment_fixed()) {
                                for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::AWAKE); i_chain++) {
                                    moment->increment_moment(MCType::AWAKE, _mc_chains.at(MCType::AWAKE).at(i_chain).at(layer2).at(sp2)(mit.row()) * _mc_chains.at(MCType::AWAKE).at(i_chain).at(layer1).at(sp1)(mit.col()) / _no_markov_chains.at(MCType::AWAKE) );
                                };
                            };
                        
                            // Iterate over asleep chains
                            for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::ASLEEP); i_chain++) {
                                moment->increment_moment(MCType::ASLEEP, _mc_chains.at(MCType::ASLEEP).at(i_chain).at(layer2).at(sp2)(mit.row()) * _mc_chains.at(MCType::ASLEEP).at(i_chain).at(layer1).at(sp1)(mit.col()) / _no_markov_chains.at(MCType::ASLEEP) );
                            };
                        };
                    };
                };
            };
        };
        
        // Reap biases
        // All biases
        int layer;
        Sptr sp;
        for (auto &bias_layer: _bias_dict) {
            layer = bias_layer.first;
            for (auto &sp_pr: bias_layer.second) {
                sp = sp_pr.first;
                moment = sp_pr.second->get_moment();
                
                // Awake phase
                
                if (!moment->get_is_awake_moment_fixed()) {
                    moment->reset_moment(MCType::AWAKE);
                    for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::AWAKE); i_chain++) {
                        moment->increment_moment(MCType::AWAKE, arma::accu(_mc_chains.at(MCType::AWAKE).at(i_chain).at(layer).at(sp)) / _no_markov_chains.at(MCType::AWAKE));
                    };
                };
                
                // Asleep phase
                
                moment->reset_moment(MCType::ASLEEP);
                for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::ASLEEP); i_chain++) {
                    moment->increment_moment(MCType::ASLEEP, arma::accu(_mc_chains.at(MCType::ASLEEP).at(i_chain).at(layer).at(sp)) / _no_markov_chains.at(MCType::ASLEEP));
                };
            };
        };
    };
    
    void Lattice::reap_moments_and_slide_centers_normal_w_centered_gradient_vec(bool slide_means, double sliding_factor) {
        
        // If centering mode: calculate the centers from the awake moment, and then slide
        if (slide_means) {
            for (auto layer=0; layer<_no_layers; layer++) {
                // Determine means
                for (auto sp: _species_possible_vec.at(layer)) {
                    // Reset
                    _c_batch_means[layer][sp].fill(arma::fill::zeros);
                    
                    // Get batch mean from all chains
                    for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::AWAKE); i_chain++) {
                        _c_batch_means[layer][sp] += _mc_chains.at(MCType::AWAKE).at(i_chain).at(layer).at(sp);
                    };
                    _c_batch_means[layer][sp] /= _no_markov_chains.at(MCType::AWAKE);
                };
                
                // Slide
                for (auto sp: _species_possible_vec.at(layer)) {
                    // Slide
                    _c_means[layer][sp] = (1.0 - sliding_factor) * _c_means.at(layer).at(sp) + sliding_factor * _c_batch_means.at(layer).at(sp);
                };
            };
        };
        
        // Reap ixns
        int layer1, layer2;
        Sptr sp1, sp2;
        arma::sp_mat::iterator mit, mit_end;
        std::shared_ptr<Moment> moment;
        arma::sp_mat* wm = nullptr;
        for (auto &o2_ixn_layer_1: _o2_ixn_dict) {
            layer1 = o2_ixn_layer_1.first;
            for (auto &sp_pr_1: o2_ixn_layer_1.second) {
                sp1 = sp_pr_1.first;
                for (auto &o2_ixn_layer_2: sp_pr_1.second) {
                    layer2 = o2_ixn_layer_2.first;
                    // Be careful not to double count
                    if (layer1 >= layer2) {
                        // skip
                        continue;
                    };
                    // Now layer1 < layer2 guaranteed
                    
                    for (auto &sp_pr_2: o2_ixn_layer_2.second) {
                        
                        sp2 = sp_pr_2.first;
                        moment = sp_pr_2.second->get_moment();
                    
                        // Awake
                        if (!moment->get_is_awake_moment_fixed()) {
                            wm = moment->get_weight_matrix_ptr(MCType::AWAKE);
                            moment->reset_weight_matrix(MCType::AWAKE);
                            for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::AWAKE); i_chain++) {
                                mit = _adj.at(layer1).at(layer2).begin();
                                mit_end = _adj.at(layer1).at(layer2).end();
                                for(; mit != mit_end; ++mit) {
                                    (*wm)(mit.row(),mit.col()) += ( _mc_chains.at(MCType::AWAKE).at(i_chain).at(layer2).at(sp2)(mit.row()) - _c_means.at(layer2).at(sp2)(mit.row()) ) * ( _mc_chains.at(MCType::AWAKE).at(i_chain).at(layer1).at(sp1)(mit.col()) - _c_means.at(layer1).at(sp1)(mit.col()) ) / _no_markov_chains.at(MCType::AWAKE);
                                };
                                // moment->increment_weight_matrix(MCType::AWAKE, _adj.at(layer1).at(layer2) % ( (_mc_chains.at(MCType::AWAKE).at(i_chain).at(layer2).at(sp2) - _c_means.at(layer2).at(sp2)) * (_mc_chains.at(MCType::AWAKE).at(i_chain).at(layer1).at(sp1).t() - _c_means.at(layer1).at(sp1).t()) ) / _no_markov_chains.at(MCType::AWAKE));
                            };
                            moment->set_moment_to_weight_matrix_sum(MCType::AWAKE);
                        };
                        
                        // Asleep
                        wm = moment->get_weight_matrix_ptr(MCType::AWAKE);
                        moment->reset_weight_matrix(MCType::ASLEEP);
                        for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::ASLEEP); i_chain++) {
                            mit = _adj.at(layer1).at(layer2).begin();
                            mit_end = _adj.at(layer1).at(layer2).end();
                            for(; mit != mit_end; ++mit) {
                                (*wm)(mit.row(),mit.col()) += ( _mc_chains.at(MCType::ASLEEP).at(i_chain).at(layer2).at(sp2)(mit.row()) - _c_means.at(layer2).at(sp2)(mit.row()) ) * ( _mc_chains.at(MCType::ASLEEP).at(i_chain).at(layer1).at(sp1)(mit.col()) - _c_means.at(layer1).at(sp1)(mit.col()) ) / _no_markov_chains.at(MCType::ASLEEP);
                            };
                            // moment->increment_weight_matrix(MCType::ASLEEP, _adj.at(layer1).at(layer2) % ( (_mc_chains.at(MCType::ASLEEP).at(i_chain).at(layer2).at(sp2) - _c_means.at(layer2).at(sp2)) * (_mc_chains.at(MCType::ASLEEP).at(i_chain).at(layer1).at(sp1).t() - _c_means.at(layer1).at(sp1).t()) ) / _no_markov_chains.at(MCType::ASLEEP));
                        };
                        moment->set_moment_to_weight_matrix_sum(MCType::ASLEEP);
                        
                        // Calculate diff
                        moment->calculate_weight_matrix_awake_minus_asleep();
                    };
                };
            };
        };

        // Reap biases
        // All biases
        int layer;
        Sptr sp;
        for (auto &bias_layer: _bias_dict) {
            layer = bias_layer.first;
            for (auto &sp_pr: bias_layer.second) {
                sp = sp_pr.first;
                moment = sp_pr.second->get_moment();
                
                // Awake phase
                
                if (!moment->get_is_awake_moment_fixed()) {
                    moment->reset_moment(MCType::AWAKE);
                    for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::AWAKE); i_chain++) {
                        moment->increment_moment(MCType::AWAKE, arma::accu(_mc_chains.at(MCType::AWAKE).at(i_chain).at(layer).at(sp)) / _no_markov_chains.at(MCType::AWAKE));
                    };
                };
                
                // Asleep phase
                
                moment->reset_moment(MCType::ASLEEP);
                for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::ASLEEP); i_chain++) {
                    moment->increment_moment(MCType::ASLEEP, arma::accu(_mc_chains.at(MCType::ASLEEP).at(i_chain).at(layer).at(sp)) / _no_markov_chains.at(MCType::ASLEEP));
                };
            };
        };
        
        // Offset bias for centering
        double offset;
        arma::vec mean;
        Iptr ixn;
        // Go through all layers
        for (auto layer=0; layer<_no_layers; layer++) {
            
            // Go through all species in this layer
            for (auto sp: _species_possible_vec.at(layer)) {
                
                // Reset offset
                offset = 0.0;
                
                // Below
                if (layer != 0) {
                    
                    // Go through all species in the layer below
                    for (auto sp_below: _species_possible_vec.at(layer-1)) {
                        
                        // Check that such an ixn exists
                        auto it1 = _o2_ixn_dict.at(layer-1).at(sp_below).find(layer);
                        if (it1 == _o2_ixn_dict.at(layer-1).at(sp_below).end()) {
                            continue;
                        };
                        auto it2 = it1->second.find(sp);
                        if (it2 == it1->second.end()) {
                            continue;
                        };
                        ixn = it2->second;
                        
                        // Calculate offset
                        offset -= arma::accu( ixn->get_moment()->get_weight_matrix_awake_minus_asleep() * _c_means.at(layer-1).at(sp_below) );
                    };
                };
                
                // Above
                if (layer != _no_layers-1) {
                    
                    // Go through all species in the layer above
                    for (auto sp_above: _species_possible_vec.at(layer+1)) {
                        
                        // Check that such an ixn exists
                        auto it1 = _o2_ixn_dict.at(layer+1).at(sp_above).find(layer);
                        if (it1 == _o2_ixn_dict.at(layer+1).at(sp_above).end()) {
                            continue;
                        };
                        auto it2 = it1->second.find(sp);
                        if (it2 == it1->second.end()) {
                            continue;
                        };
                        ixn = it2->second;
                        
                        // Calculate offset
                        offset -= arma::accu( ixn->get_moment()->get_weight_matrix_awake_minus_asleep().t() * _c_means.at(layer+1).at(sp_above) );
                    };
                };
                
                // Adjust bias in this layer
                _bias_dict.at(layer).at(sp)->get_moment()->set_moment_diff_awake_minus_asleep_offset(offset);
            };
        };

    };
    
    void Lattice::reap_moments_and_slide_centers_normal_w_centered_gradient_pt(bool slide_means, double sliding_factor) {
        
        // If centering mode: calculate the centers from the awake moment, and then slide
        if (slide_means) {
            int no_units;
            for (auto layer=0; layer<_no_layers; layer++) {
                no_units = get_no_units_in_layer(layer);
                
                // Determine means
                for (auto sp: _species_possible_vec.at(layer)) {
                    // Reset
                    _cpt_batch_means[layer][sp] = 0.0;
                    
                    // Get batch mean from all chains
                    for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::AWAKE); i_chain++) {
                        _cpt_batch_means[layer][sp] += arma::accu(_mc_chains.at(MCType::AWAKE).at(i_chain).at(layer).at(sp)) / no_units;
                    };
                    _cpt_batch_means[layer][sp] /= _no_markov_chains.at(MCType::AWAKE);
                };
                
                // Slide
                for (auto sp: _species_possible_vec.at(layer)) {
                    // Slide
                    _cpt_means[layer][sp] = (1.0 - sliding_factor) * _cpt_means.at(layer).at(sp) + sliding_factor * _cpt_batch_means.at(layer).at(sp);
                };
            };
        };
        
        // Reap ixns
        int layer1, layer2;
        Sptr sp1, sp2;
        arma::sp_mat::iterator mit, mit_end;
        std::shared_ptr<Moment> moment;
        for (auto &o2_ixn_layer_1: _o2_ixn_dict) {
            layer1 = o2_ixn_layer_1.first;
            for (auto &sp_pr_1: o2_ixn_layer_1.second) {
                sp1 = sp_pr_1.first;
                for (auto &o2_ixn_layer_2: sp_pr_1.second) {
                    layer2 = o2_ixn_layer_2.first;
                    // Be careful not to double count
                    if (layer1 >= layer2) {
                        // skip
                        continue;
                    };
                    // Now layer1 < layer2 guaranteed
                    
                    for (auto &sp_pr_2: o2_ixn_layer_2.second) {
                        
                        sp2 = sp_pr_2.first;
                        moment = sp_pr_2.second->get_moment();
                        
                        // Awake
                        if (!moment->get_is_awake_moment_fixed()) {
                            moment->reset_moment(MCType::AWAKE);
                            for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::AWAKE); i_chain++) {
                                mit = _adj.at(layer1).at(layer2).begin();
                                mit_end = _adj.at(layer1).at(layer2).end();
                                for(; mit != mit_end; ++mit) {
                                    moment->increment_moment(MCType::AWAKE, ( _mc_chains.at(MCType::AWAKE).at(i_chain).at(layer2).at(sp2)(mit.row()) - _cpt_means.at(layer2).at(sp2) ) * ( _mc_chains.at(MCType::AWAKE).at(i_chain).at(layer1).at(sp1)(mit.col()) - _cpt_means.at(layer1).at(sp1) ) / _no_markov_chains.at(MCType::AWAKE) );
                                };
                            };
                        };
                        
                        // Asleep phase
                        moment->reset_moment(MCType::ASLEEP);
                        for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::ASLEEP); i_chain++) {
                            mit = _adj.at(layer1).at(layer2).begin();
                            mit_end = _adj.at(layer1).at(layer2).end();
                            for(; mit != mit_end; ++mit) {
                                moment->increment_moment(MCType::ASLEEP, ( _mc_chains.at(MCType::ASLEEP).at(i_chain).at(layer2).at(sp2)(mit.row()) - _cpt_means.at(layer2).at(sp2) ) * ( _mc_chains.at(MCType::ASLEEP).at(i_chain).at(layer1).at(sp1)(mit.col()) - _cpt_means.at(layer1).at(sp1) ) / _no_markov_chains.at(MCType::ASLEEP) );
                            };
                        };
                    };
                };
            };
        };
        
        // Reap biases
        // All biases
        int layer;
        Sptr sp;
        for (auto &bias_layer: _bias_dict) {
            layer = bias_layer.first;
            for (auto &sp_pr: bias_layer.second) {
                sp = sp_pr.first;
                moment = sp_pr.second->get_moment();
                
                // Awake phase
                
                if (!moment->get_is_awake_moment_fixed()) {
                    moment->reset_moment(MCType::AWAKE);
                    for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::AWAKE); i_chain++) {
                        moment->increment_moment(MCType::AWAKE, arma::accu(_mc_chains.at(MCType::AWAKE).at(i_chain).at(layer).at(sp)) / _no_markov_chains.at(MCType::AWAKE));
                    };
                };
                
                // Asleep phase
                
                moment->reset_moment(MCType::ASLEEP);
                for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::ASLEEP); i_chain++) {
                    moment->increment_moment(MCType::ASLEEP, arma::accu(_mc_chains.at(MCType::ASLEEP).at(i_chain).at(layer).at(sp)) / _no_markov_chains.at(MCType::ASLEEP));
                };
            };
        };
        
        // Offset bias for centering
        double offset;
        arma::vec mean;
        Iptr ixn;
        // Go through all layers
        for (auto layer=0; layer<_no_layers; layer++) {
            
            // Go through all species in this layer
            for (auto sp: _species_possible_vec.at(layer)) {
                
                // Reset offset
                offset = 0.0;
                
                // Below
                if (layer != 0) {
                    
                    // Go through all species in the layer below
                    for (auto sp_below: _species_possible_vec.at(layer-1)) {
                        
                        // Check that such an ixn exists
                        auto it1 = _o2_ixn_dict.at(layer-1).at(sp_below).find(layer);
                        if (it1 == _o2_ixn_dict.at(layer-1).at(sp_below).end()) {
                            continue;
                        };
                        auto it2 = it1->second.find(sp);
                        if (it2 == it1->second.end()) {
                            continue;
                        };
                        ixn = it2->second;
                        
                        offset -= 8.0 * ixn->get_moment()->get_moment_diff_awake_minus_asleep() * _cpt_means.at(layer-1).at(sp_below);
                    };
                };

                // Above
                if (layer != _no_layers-1) {
                    
                    // Go through all species in the layer above
                    for (auto sp_above: _species_possible_vec.at(layer+1)) {
                        
                        // Check that such an ixn exists
                        auto it1 = _o2_ixn_dict.at(layer+1).at(sp_above).find(layer);
                        if (it1 == _o2_ixn_dict.at(layer+1).at(sp_above).end()) {
                            continue;
                        };
                        auto it2 = it1->second.find(sp);
                        if (it2 == it1->second.end()) {
                            continue;
                        };
                        ixn = it2->second;
                        
                        // Calculate offset
                        offset -= 8.0 * ixn->get_moment()->get_moment_diff_awake_minus_asleep() * _cpt_means.at(layer+1).at(sp_above);
                    };
                };
                
                // Adjust bias in this layer
                _bias_dict.at(layer).at(sp)->get_moment()->set_moment_diff_awake_minus_asleep_offset(offset);
            };
        };

    };
    
    void Lattice::reap_moments_and_slide_centers_centered_pt(bool slide_means, double sliding_factor) {

        // If centered pt M mode: determine means, adjust bias, THEN slide!
        if (slide_means) {
            int no_units;
            for (auto layer=0; layer<_no_layers; layer++) {
                
                no_units = get_no_units_in_layer(layer);
                
                // Determine means
                for (auto sp: _species_possible_vec.at(layer)) {
                    // Reset
                    _cpt_batch_means[layer][sp] = 0.0;
                    
                    // Get batch mean from all chains
                    for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::AWAKE); i_chain++) {
                        _cpt_batch_means[layer][sp] += arma::accu(_mc_chains.at(MCType::AWAKE).at(i_chain).at(layer).at(sp)) / no_units;
                    };
                    _cpt_batch_means[layer][sp] /= _no_markov_chains.at(MCType::AWAKE);
                };
            };
            double val;
            for (auto layer=0; layer<_no_layers; layer++) {
                for (auto sp: _species_possible_vec.at(layer)) {
                    val = _bias_dict.at(layer).at(sp)->get_val();
                    
                    // Layer below
                    for (auto sp_below: _species_possible_vec.at(layer-1)) {
                        val += 8.0 * sliding_factor * _o2_ixn_dict.at(layer-1).at(sp_below).at(layer).at(sp)->get_val() * ( _cpt_batch_means.at(layer-1).at(sp_below) - _cpt_means.at(layer-1).at(sp_below) );
                    };
                
                    // Layer above
                    for (auto sp_above: _species_possible_vec.at(layer+1)) {
                        val += 8.0 * sliding_factor * _o2_ixn_dict.at(layer+1).at(sp_above).at(layer).at(sp)->get_val() * ( _cpt_batch_means.at(layer+1).at(sp_above) - _cpt_means.at(layer+1).at(sp_above) );
                    };
                
                    // Set val
                    _bias_dict.at(layer).at(sp)->set_val(val);
                };
            };
            for (auto layer=0; layer<_no_layers; layer++) {
                for (auto sp: _species_possible_vec.at(layer)) {
                    // Slide
                    _cpt_means[layer][sp] = (1.0 - sliding_factor) * _cpt_means.at(layer).at(sp) + sliding_factor * _cpt_batch_means.at(layer).at(sp);
                };
            };
        };
        
        // Reap ixns
        int layer1, layer2;
        Sptr sp1, sp2;
        arma::sp_mat::iterator mit, mit_end;
        std::shared_ptr<Moment> moment;
        for (auto &o2_ixn_layer_1: _o2_ixn_dict) {
            layer1 = o2_ixn_layer_1.first;
            for (auto &sp_pr_1: o2_ixn_layer_1.second) {
                sp1 = sp_pr_1.first;
                for (auto &o2_ixn_layer_2: sp_pr_1.second) {
                    layer2 = o2_ixn_layer_2.first;
                    // Be careful not to double count
                    if (layer1 >= layer2) {
                        // skip
                        continue;
                    };
                    // Now layer1 < layer2 guaranteed
                    
                    for (auto &sp_pr_2: o2_ixn_layer_2.second) {
                        
                        sp2 = sp_pr_2.first;
                        moment = sp_pr_2.second->get_moment();
                        
                        // Awake
                        if (!moment->get_is_awake_moment_fixed()) {
                            moment->reset_moment(MCType::AWAKE);
                            for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::AWAKE); i_chain++) {
                                mit = _adj.at(layer1).at(layer2).begin();
                                mit_end = _adj.at(layer1).at(layer2).end();
                                for(; mit != mit_end; ++mit) {
                                    moment->increment_moment(MCType::AWAKE, ( _mc_chains.at(MCType::AWAKE).at(i_chain).at(layer2).at(sp2)(mit.row()) - _cpt_means.at(layer2).at(sp2) ) * ( _mc_chains.at(MCType::AWAKE).at(i_chain).at(layer1).at(sp1)(mit.col()) - _cpt_means.at(layer1).at(sp1) ) / _no_markov_chains.at(MCType::AWAKE) );
                                };
                            };
                        };
                        
                        // Asleep phase
                        moment->reset_moment(MCType::ASLEEP);
                        for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::ASLEEP); i_chain++) {
                            mit = _adj.at(layer1).at(layer2).begin();
                            mit_end = _adj.at(layer1).at(layer2).end();
                            for(; mit != mit_end; ++mit) {
                                moment->increment_moment(MCType::ASLEEP, ( _mc_chains.at(MCType::ASLEEP).at(i_chain).at(layer2).at(sp2)(mit.row()) - _cpt_means.at(layer2).at(sp2) ) * ( _mc_chains.at(MCType::ASLEEP).at(i_chain).at(layer1).at(sp1)(mit.col()) - _cpt_means.at(layer1).at(sp1) ) / _no_markov_chains.at(MCType::ASLEEP) );
                            };
                        };
                    };
                };
            };
        };
        
        // Reap biases
        // All biases
        int layer;
        Sptr sp;
        for (auto &bias_layer: _bias_dict) {
            layer = bias_layer.first;
            for (auto &sp_pr: bias_layer.second) {
                sp = sp_pr.first;
                moment = sp_pr.second->get_moment();
                
                // Awake phase
                
                if (!moment->get_is_awake_moment_fixed()) {
                    moment->reset_moment(MCType::AWAKE);
                    for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::AWAKE); i_chain++) {
                        moment->increment_moment(MCType::AWAKE, arma::accu(_mc_chains.at(MCType::AWAKE).at(i_chain).at(layer).at(sp)) / _no_markov_chains.at(MCType::AWAKE));
                    };
                };
                
                // Asleep phase
                
                moment->reset_moment(MCType::ASLEEP);
                for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::ASLEEP); i_chain++) {
                    moment->increment_moment(MCType::ASLEEP, arma::accu(_mc_chains.at(MCType::ASLEEP).at(i_chain).at(layer).at(sp)) / _no_markov_chains.at(MCType::ASLEEP));
                };
            };
        };
    };
    
    // ***************
    // MARK: - Write out centers
    // ***************
    
    void Lattice::read_centers_from_file(int layer, std::string fname) {
        // Open
        std::ifstream f;
        f.open(fname);
        if (!f.is_open()) { // make sure we found it
            std::cerr << ">>> Error: Lattice::read_centers_from_file <<< could not find file: " << fname << std::endl;
            exit(EXIT_FAILURE);
        };
        
        std::map<std::string,std::map<int,double>> to_write;
        
        std::string x="",y="",z="";
        std::string sp="";
        std::string center="";
        std::string line;
        std::istringstream iss;
        int s;
        
        int no_species = _species_possible_vec.at(layer).size();
        
        if (_no_dims == 1) {
            while (getline(f,line)) {
                if (line == "") { continue; };
                iss = std::istringstream(line);
                iss >> x;
                s = _look_up_unit(layer,atoi(x.c_str()));
                for (auto i=0; i<no_species; i++) {
                    iss >> sp >> center;
                    to_write[sp][s] = atof(center.c_str());;
                    sp=""; center="";
                };
                x="";
            };
        } else if (_no_dims == 2) {
            while (getline(f,line)) {
                if (line == "") { continue; };
                iss = std::istringstream(line);
                iss >> x >> y;
                s = _look_up_unit(layer,atoi(x.c_str()),atoi(y.c_str()));
                for (auto i=0; i<no_species; i++) {
                    iss >> sp >> center;
                    to_write[sp][s] = atof(center.c_str());;
                    sp=""; center="";
                };
                x=""; y="";
            };
        }else if (_no_dims == 3) {
            while (getline(f,line)) {
                if (line == "") { continue; };
                iss = std::istringstream(line);
                iss >> x >> y >> z;
                s = _look_up_unit(layer,atoi(x.c_str()),atoi(y.c_str()),atoi(z.c_str()));
                for (auto i=0; i<no_species; i++) {
                    iss >> sp >> center;
                    to_write[sp][s] = atof(center.c_str());;
                    sp=""; center="";
                };
                x=""; y=""; z="";
            };
        };
        
        // Write
        for (auto sp_pr: to_write) {
            auto species = _species_possible_map[layer][sp_pr.first];
            for (auto idx_pr: sp_pr.second) {
                _c_means[layer][species](idx_pr.first) = idx_pr.second;
            };
        };
        
        // Close!!!
        f.close();
    };
    
    void Lattice::write_centers_to_file(int layer, std::string fname) const {
        std::ofstream f;
        f.open (fname);
        if (!f.is_open()) { // make sure we found it
            std::cerr << ">>> Lattice::write_centers_to_file <<< Error: could not open file: " << fname << " for writing" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        // Go through all units
        int no_units = get_no_units_in_layer(layer);

        // Go through all units
        std::vector<int> pos;
        for (auto unit=0; unit<no_units; unit++) {
            
            // Go through species possible
            pos = _look_up_pos(layer,unit);
            for (auto const &x: pos) {
                f << x << " ";
            };
            
            // Write species and mean
            for (auto sp: _species_possible_vec.at(layer)) {
                f << sp->get_name() << " " << _c_means.at(layer).at(sp)(unit) << " ";
            };
            f << "\n";
         };
        
        // Close!!!
        f.close();
    };
    
    void Lattice::read_center_pts_from_file(std::string fname) {
        // Open
        std::ifstream f;
        f.open(fname);
        if (!f.is_open()) { // make sure we found it
            std::cerr << ">>> Error: Lattice::read_center_pt_from_file <<< could not find file: " << fname << std::endl;
            exit(EXIT_FAILURE);
        };
        
        std::string sp="", center="", layer_str="";
        std::string line;
        std::istringstream iss;
        Sptr species;
        int layer;

        while (getline(f,line)) {
            if (line == "") { continue; };
            iss = std::istringstream(line);
            iss >> layer_str >> sp >> center;
            if (sp != "") {
                layer = atoi(layer_str.c_str());
                species = _species_possible_map.at(layer).at(sp);
                _cpt_means[layer][species] = atof(center.c_str());
                // std::cout << "Lattice::read_center_pt_from_file: Set layer: " << layer << " species: " << species->get_name() << " to center: " << _cpt_means.at(layer).at(species) << std::endl;
            };
            sp=""; center=""; layer_str="";
        };
        
        // Close!!!
        f.close();
    };
    
    void Lattice::write_center_pts_to_file(std::string fname) const {
        std::ofstream f;
        f.open (fname);
        if (!f.is_open()) { // make sure we found it
            std::cerr << ">>> Lattice::write_center_pt_to_file <<< Error: could not open file: " << fname << " for writing" << std::endl;
            exit(EXIT_FAILURE);
        };

        for (auto layer=0; layer<_no_layers; layer++) {
            for (auto pr: _cpt_means.at(layer)) {
                f << layer << " " << pr.first->get_name() << " " << pr.second << "\n";
            };
        };
        
        // Close!!!
        f.close();
    };
    
    // ***************
    // MARK: - Set centers
    // ***************
    
    double Lattice::get_center_pt_for_species_in_layer(int layer, Sptr species) const {
        return _cpt_means.at(layer).at(species);
    };
    void Lattice::set_center_pt_for_species_in_layer(int layer, Sptr species, double center) {
        _cpt_means[layer][species] = center;
    };

    arma::vec Lattice::get_center_vec_for_species_in_layer(int layer, Sptr species) const {
        return _c_means.at(layer).at(species);
    };
    void Lattice::set_center_vec_for_species_in_layer(int layer, Sptr species, arma::vec center) {
        _c_means[layer][species] = center;
    };
    void Lattice::set_center_vec_for_species_in_layer(int layer, Sptr species, double center) {
        int no_units = get_no_units_in_layer(layer);
        _c_means[layer][species] = center * arma::vec(no_units,arma::fill::ones);
    };
    
    // ***************
    // MARK: - Add ixn to all ixns vec
    // ***************
    
    // Add ixn to all ixns vec
    void Lattice::_add_to_all_ixns_vec(Iptr ixn) {
    
        auto it = std::find(_all_ixns.begin(),_all_ixns.end(),ixn);
        if (it == _all_ixns.end()) {
            _all_ixns.push_back(ixn);
        } else {
            /*
            std::cerr << ">>> Lattice::_add_to_all_ixns_vec <<< reusing ixns currently not supported because it is not treated correctly in the reap function!" << std::endl;
            exit(EXIT_FAILURE);
             */
        };
    };
    
    // ***************
    // MARK: - Wake/sleep loop
    // ***************
    
    void Lattice::wake_sleep_loop_bm_pcd(int i_opt_step, int no_mean_field_updates, int no_gibbs_sampling_steps, std::vector<FName> &fnames, OptionsWakeSleep_BM_PCD options) {
        
        // AWAKE PHASE
        
        clock_t t0 = clock();
        
        // Read in the batch
        for (int i_chain=0; i_chain<_no_markov_chains[MCType::AWAKE]; i_chain++)
        {
            read_layer_from_file(MCType::AWAKE, i_chain, 0, fnames[i_chain].name, fnames[i_chain].binary);
        };
        
        clock_t t1 = clock();
        
        // Option (1): init MF with random hidden layers with prob units
        /*
         for (int i_chain=0; i_chain<_no_markov_chains[MCType::AWAKE]; i_chain++) {
         _latt->set_random_all_hidden_units(MCType::AWAKE, i_chain, false);
         };
         */
        // Option (2): upward pass with 2x weights (DBM) to activate probabilitsic units
        // (faster to converge!!!)
        activate_upward_pass_with_2x_weights_1x_bias(MCType::AWAKE, false);
        
        clock_t t2 = clock();
        
        // Variational inference
        if (!options.gibbs_sample_awake_phase) {
            
            for (auto i=0; i<no_mean_field_updates; i++) {
                mean_field_hiddens_step();
            };
            // gibbs_sampling_step_awake(true);
            
        } else {
            
            for (auto i=0; i<no_mean_field_updates; i++) {
                gibbs_sampling_step_awake(options.gibbs_sample_awake_phase_hidden_binary);
            };
        };
        
        clock_t t3 = clock();
        
        // Write out the lattices
        if (options.write_after_awake) {
            for (auto i_chain=0; i_chain<_no_markov_chains[MCType::AWAKE]; i_chain++) {
                for (auto layer=0; layer<_no_layers; layer++) {
                    bool write_binary = true;
                    if (layer != 0) {
                        write_binary = false;
                    };
                    write_layer_to_file(MCType::AWAKE, i_chain, layer, options.write_after_awake_dir+"/"+pad_str(i_chain,2)+"_"+pad_str(layer,2)+".txt", write_binary);
                };
            };
        };
        
        // ASLEEP PHASE - PERSISTENT_CD
        
        // Run CD sampling
        
        // Sample vis, hidden
        for (int i_sampling_step=0; i_sampling_step<no_gibbs_sampling_steps-1; i_sampling_step++)
        {
            gibbs_sampling_step(options.is_asleep_visible_binary, options.is_asleep_hidden_binary);
        };
        // Final step
        if (options.is_asleep_visible_binary_final && options.is_asleep_hidden_binary_final) {
            // All binary
            gibbs_sampling_step(options.is_asleep_visible_binary_final, options.is_asleep_hidden_binary_final);
        } else {
            // Parallel for non-binary options
            gibbs_sampling_step_parallel(options.is_asleep_visible_binary_final, options.is_asleep_hidden_binary_final);
        };
        
        clock_t t4 = clock();
        
        // Write out the lattices
        if (options.write_after_asleep) {
            for (auto i_chain=0; i_chain<_no_markov_chains[MCType::ASLEEP]; i_chain++) {
                for (auto layer=0; layer<_no_layers; layer++) {
                    bool write_binary = true;
                    if (layer != 0) {
                        write_binary = false;
                    };
                    write_layer_to_file(MCType::ASLEEP, i_chain, layer, options.write_after_asleep_dir+"/"+pad_str(i_chain,2)+"_"+pad_str(layer,2)+".txt", write_binary);
                };
            };
        };
        
        double dt1 = (t1-t0)  / (double) CLOCKS_PER_SEC;
        double dt2 = (t2-t1)  / (double) CLOCKS_PER_SEC;
        double dt3 = (t3-t2)  / (double) CLOCKS_PER_SEC;
        double dt4 = (t4-t3)  / (double) CLOCKS_PER_SEC;
        double dt_tot = dt1 + dt2 + dt3 + dt4;
        if (options.verbose_timing) {
            std::cout << "[time " << dt_tot << "] [read " << dt1/dt_tot << "] [up " << dt2/dt_tot << "] [mf " << dt3/dt_tot << "] [gibbs " << dt4/dt_tot << "]" << std::endl;
        };
    };
    
    void Lattice::wake_sleep_loop_rbm_cd(int i_opt_step, int no_cd_steps, std::vector<FName> &fnames, OptionsWakeSleep_RBM_CD options) {
        
        // AWAKE PHASE
        
        clock_t t0 = clock();
        
        // Read in the batch
        for (int i_chain=0; i_chain<_no_markov_chains[MCType::AWAKE]; i_chain++)
        {
            read_layer_from_file(MCType::AWAKE, i_chain, 0, fnames[i_chain].name, fnames[i_chain].binary);
            // read_layer_from_file(MCType::ASLEEP, i_chain, 0, fnames[i_chain].name, fnames[i_chain].binary);
        };
        _mc_chains[MCType::ASLEEP] = _mc_chains.at(MCType::AWAKE);
        
        clock_t t1 = clock();
        
        // AWAKE PHASE
        
        // Activate hidden layer; use probs!
        activate_layer_calculate_from_below(MCType::AWAKE, 1);
        activate_layer_convert_to_probs(MCType::AWAKE, 1, false);
        activate_layer_committ(MCType::AWAKE, 1);

        clock_t t2 = clock();
        
        // ASLEEP PHASE - PERSISTENT_CD
    
        for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::ASLEEP); i_chain++) {
            set_random_all_units_in_layer(MCType::ASLEEP, i_chain, 0, true);
        };
        
        // Run CD sampling
    
        // Activate the hiddens with binary
        activate_layer_calculate_from_below(MCType::ASLEEP, 1);
        activate_layer_convert_to_probs(MCType::ASLEEP, 1, true);
        activate_layer_committ(MCType::ASLEEP, 1);
        
        // Sample vis, hidden repeadedly
        for (int i_sampling_step=0; i_sampling_step<no_cd_steps-1; i_sampling_step++)
        {
            activate_layer_calculate_from_above(MCType::ASLEEP, 0);
            activate_layer_convert_to_probs(MCType::ASLEEP, 0, true);
            activate_layer_committ(MCType::ASLEEP, 0);

            activate_layer_calculate_from_below(MCType::ASLEEP, 1);
            activate_layer_convert_to_probs(MCType::ASLEEP, 1, true);
            activate_layer_committ(MCType::ASLEEP, 1);
        };
        // Final step: use probs for hidden layer
        activate_layer_calculate_from_above(MCType::ASLEEP, 0);
        activate_layer_convert_to_probs(MCType::ASLEEP, 0, true);
        activate_layer_committ(MCType::ASLEEP, 0);
        
        activate_layer_calculate_from_below(MCType::ASLEEP, 1);
        activate_layer_convert_to_probs(MCType::ASLEEP, 1, false);
        activate_layer_committ(MCType::ASLEEP, 1);
        
        clock_t t3 = clock();

        if (options.verbose_timing) {
            double dt1 = (t1-t0)  / (double) CLOCKS_PER_SEC;
            double dt2 = (t2-t1)  / (double) CLOCKS_PER_SEC;
            double dt3 = (t3-t2)  / (double) CLOCKS_PER_SEC;
            double dt_tot = dt1 + dt2 + dt3;
            std::cout << "[time " << dt_tot << "] [read " << dt1/dt_tot << "] [awake " << dt2/dt_tot << "] [asleep " << dt3/dt_tot << "]" << std::endl;
        };
    };
    
    // ***************
    // MARK: - Counts
    // ***************
    
    double Lattice::get_count_1d(MCType chain, int i_chain, Sptr sp) const {
        if (_no_dims != 1) {
            std::cerr << ">>> Lattice::get_count_1d <<< Error: only for 1D lattices" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        return arma::accu(_mc_chains.at(chain).at(i_chain).at(0).at(sp));
    };
    double Lattice::get_count_1d(MCType chain, int i_chain, Sptr sp1, Sptr sp2) const {
        if (_no_dims != 1) {
            std::cerr << ">>> Lattice::get_count_1d <<< Error: only for 1D lattices" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        double count=0.;
        for (auto x=0; x<_box_length-1; x++) {
            count += _mc_chains.at(chain).at(i_chain).at(0).at(sp1).at(x) * _mc_chains.at(chain).at(i_chain).at(0).at(sp2).at(x+1);
        };
        if (sp1 != sp2) {
            for (auto x=0; x<_box_length-1; x++) {
                count += _mc_chains.at(chain).at(i_chain).at(0).at(sp2).at(x) * _mc_chains.at(chain).at(i_chain).at(0).at(sp1).at(x+1);
            };
        };
        return count;
    };
    double Lattice::get_count_1d(MCType chain, int i_chain, Sptr sp1, Sptr sp2, Sptr sp3) const {
        if (_no_dims != 1) {
            std::cerr << ">>> Lattice::get_count_1d <<< Error: only for 1D lattices" << std::endl;
            exit(EXIT_FAILURE);
        };

        double count=0.;
        for (auto x=0; x<_box_length-2; x++) {
            count += _mc_chains.at(chain).at(i_chain).at(0).at(sp1).at(x) * _mc_chains.at(chain).at(i_chain).at(0).at(sp2).at(x+1) * _mc_chains.at(chain).at(i_chain).at(0).at(sp3).at(x+2);
        };
        if (!((sp1 == sp2) && (sp2 == sp3))) {
            for (auto x=0; x<_box_length-2; x++) {
                count += _mc_chains.at(chain).at(i_chain).at(0).at(sp3).at(x) * _mc_chains.at(chain).at(i_chain).at(0).at(sp2).at(x+1) * _mc_chains.at(chain).at(i_chain).at(0).at(sp1).at(x+2);
            };
        };
        return count;
    };
    double Lattice::get_count_1d(MCType chain, int i_chain, Sptr sp1, Sptr sp2, Sptr sp3, Sptr sp4) const {
        if (_no_dims != 1) {
            std::cerr << ">>> Lattice::get_count_1d <<< Error: only for 1D lattices" << std::endl;
            exit(EXIT_FAILURE);
        };

        double count=0.;
        for (auto x=0; x<_box_length-3; x++) {
            count += _mc_chains.at(chain).at(i_chain).at(0).at(sp1).at(x) * _mc_chains.at(chain).at(i_chain).at(0).at(sp2).at(x+1) * _mc_chains.at(chain).at(i_chain).at(0).at(sp3).at(x+2) * _mc_chains.at(chain).at(i_chain).at(0).at(sp4).at(x+3);
        };
        if (!((sp1 == sp2) && (sp2 == sp3) && (sp3 == sp4))) {
            for (auto x=0; x<_box_length-3; x++) {
                count += _mc_chains.at(chain).at(i_chain).at(0).at(sp4).at(x) * _mc_chains.at(chain).at(i_chain).at(0).at(sp3).at(x+1) * _mc_chains.at(chain).at(i_chain).at(0).at(sp2).at(x+2) * _mc_chains.at(chain).at(i_chain).at(0).at(sp1).at(x+3);
            };
        };
        return count;
    };

};
