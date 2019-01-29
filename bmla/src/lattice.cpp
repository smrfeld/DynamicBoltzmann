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

	Lattice::Lattice(int dim, int box_length, std::vector<Sptr> species_visible)
	{
		if (dim != 1 && dim != 2 && dim != 3) {
			std::cerr << "ERROR: only dimensions 1,2,3 are supported for Lattice." << std::endl;
			exit(EXIT_FAILURE);
		};
		_no_dims = dim;
		_box_length = box_length;
        _no_markov_chains = 1;
        _no_markov_chains_asleep = 0;
        _no_layers = 1;
        _i_markov_chain = 0;
        
		// Visible layer
        int i_chain = 0;
        int i_layer = 0;
		if (dim == 1) {
            for (auto sp: species_visible) {
                _latt[i_chain][i_layer][sp] = arma::vec(_box_length,arma::fill::zeros);
                _latt_act[i_chain][i_layer][sp] = arma::vec(_box_length,arma::fill::zeros);
            };
            for (auto x=1; x<=_box_length; x++) {
                _lookup_1[0][x] = x-1;
                _rlookup[0][x-1] = std::vector<int>({x});
            };
        } else if (dim == 2) {
            for (auto sp: species_visible) {
                _latt[i_chain][i_layer][sp] = arma::vec(pow(_box_length,2),arma::fill::zeros);
                _latt_act[i_chain][i_layer][sp] = arma::vec(pow(_box_length,2),arma::fill::zeros);
            };
            int ctr=0;
            for (auto x=1; x<=_box_length; x++) {
                for (auto y=1; y<=_box_length; y++) {
                    _lookup_2[0][x][y] = ctr;
                    _rlookup[0][ctr++] = std::vector<int>({x,y});
                };
            };
		} else if (dim == 3) {
            for (auto sp: species_visible) {
                _latt[i_chain][i_layer][sp] = arma::vec(pow(_box_length,3),arma::fill::zeros);
                _latt_act[i_chain][i_layer][sp] = arma::vec(pow(_box_length,3),arma::fill::zeros);
            };
            int ctr=0;
            for (auto x=1; x<=_box_length; x++) {
                for (auto y=1; y<=_box_length; y++) {
                    for (auto z=1; z<=_box_length; z++) {
                        _lookup_3[0][x][y][z] = ctr;
                        _rlookup[0][ctr++] = std::vector<int>({x,y,z});
                    };
                };
            };
		};
        
        // No free idxs in first layer
        _free_idxs[0] = pow(_box_length,_no_dims);
        
        // Possible species
        for (auto sp: species_visible) {
            _species_possible[0][sp->get_name()] = sp;
        };
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
        _no_markov_chains_asleep = other._no_markov_chains_asleep;
        _no_layers = other._no_layers;
        _i_markov_chain = other._i_markov_chain;
        
        _latt = other._latt;
        _lookup_1 = other._lookup_1;
        _lookup_2 = other._lookup_2;
        _lookup_3 = other._lookup_3;
        _rlookup = other._rlookup;
        _free_idxs = other._free_idxs;
        _adj = other._adj;
        
        _all_ixns = other._all_ixns;
        _bias_dicts = other._bias_dicts;
        _o2_ixn_dicts = other._o2_ixn_dicts;
        _o2_mults = other._o2_mults;
        
        _species_possible = other._species_possible;
    };
	void Lattice::_move(Lattice& other) {
        _no_dims = other._no_dims;
        _box_length = other._box_length;
        _no_markov_chains = other._no_markov_chains;
        _no_markov_chains_asleep = other._no_markov_chains_asleep;
        _no_layers = other._no_layers;
        _i_markov_chain = other._i_markov_chain;
        
        _latt = other._latt;
        _lookup_1 = other._lookup_1;
        _lookup_2 = other._lookup_2;
        _lookup_3 = other._lookup_3;
        _rlookup = other._rlookup;
        _free_idxs = other._free_idxs;
        _adj = other._adj;
        
        _all_ixns = other._all_ixns;
        _bias_dicts = other._bias_dicts;
        _o2_ixn_dicts = other._o2_ixn_dicts;
        _o2_mults = other._o2_mults;
        
        _species_possible = other._species_possible;

		// Reset other
		other._no_dims = 0;
		other._box_length = 0;
        other._no_markov_chains = 0;
        other._no_markov_chains_asleep = 0;
        other._no_layers = 0;
        other._i_markov_chain = 0;
        
        other._latt.clear();
		other._lookup_1.clear();
		other._lookup_2.clear();
		other._lookup_3.clear();
        other._rlookup.clear();
        other._free_idxs.clear();
        other._adj.clear();
        
        other._all_ixns.clear();
        other._bias_dicts.clear();
        other._o2_ixn_dicts.clear();
        other._o2_mults.clear();
        
        other._species_possible.clear();
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
        auto it = _latt.at(0).at(layer).begin();
        if (it != _latt.at(0).at(layer).end()) {
            return it->second.size();
        } else {
            return 0;
        };
    };
    
    /********************
     Markov chains
     ********************/
    
    int Lattice::get_no_markov_chains_asleep() const {
        return _no_markov_chains_asleep;
    };
    void Lattice::set_no_markov_chains_asleep(int no_markov_chains_asleep) {
        _no_markov_chains_asleep = no_markov_chains_asleep;
        _no_markov_chains = _no_markov_chains_asleep + 1;
        
        if (_latt.size() < _no_markov_chains) {
            for (auto i_chain=_latt.size(); i_chain < _no_markov_chains; i_chain++) {
                _latt[i_chain] = _latt[i_chain-1];
                _latt_act[i_chain] = _latt_act[i_chain-1];
            };
        };
        if (_latt.size() > _no_markov_chains) {
            for (auto i_chain=_latt.size(); i_chain > _no_markov_chains; i_chain--) {
                _latt.erase(i_chain);
                _latt_act.erase(i_chain);
            };
        };
    };
    void Lattice::switch_to_markov_chain_asleep(int i_markov_chain_asleep) {
        if (i_markov_chain_asleep > _no_markov_chains_asleep) {
            std::cerr << ">>> Lattice::switch_to_markov_chain_asleep <<< there are only: " << _no_markov_chains_asleep << " asleep markov chains, but you tried to switch to: " << i_markov_chain_asleep << std::endl;
            exit(EXIT_FAILURE);
        };
        _i_markov_chain = i_markov_chain_asleep+1;
    };
    void Lattice::switch_to_markov_chain_awake() {
        _i_markov_chain = 0;
    };
    
    /********************
     Add a layer
     ********************/
    
    void Lattice::add_layer(int layer, int no_units, std::vector<Sptr> species) {
        if (layer != _no_layers) {
            std::cerr << ">>> Lattice::add_layer <<< error: next layer must be: " << _no_layers << " not: " << layer << std::endl;
            exit(EXIT_FAILURE);
        };
        
        // Increment no layers
        _no_layers += 1;
        
        // Add empties
        for (auto i_chain=0; i_chain<_no_markov_chains; i_chain++) {
            for (auto sp: species) {
                _latt[i_chain][layer][sp] = arma::vec(no_units,arma::fill::zeros);
                _latt_act[i_chain][layer][sp] = arma::vec(no_units,arma::fill::zeros);
            };
        };
        
        // Free idxs
        _free_idxs[layer] = 0;
        
        // Add species possible
        for (auto sp: species) {
            _species_possible[layer][sp->get_name()] = sp;
        };
        
        // Add adjacency matrix
        int size_below = get_no_units_in_layer(layer-1);
        _adj[layer-1] = arma::mat(size_below,no_units,arma::fill::zeros);
    };
    
    void Lattice::set_pos_of_hidden_unit(int layer, int x) {
        if (_free_idxs[layer] >= get_no_units_in_layer(layer)) {
            std::cerr << ">>> Lattice::set_pos_of_hidden_unit <<< no more free units in layer: " << layer << " of size: " << get_no_units_in_layer(layer) << std::endl;
            exit(EXIT_FAILURE);
        };
        
        int new_idx = _free_idxs[layer];
        _lookup_1[layer][x] = new_idx;
        _rlookup[layer][new_idx] = std::vector<int>({x});
        _free_idxs[layer] += 1;
    };
    void Lattice::set_pos_of_hidden_unit(int layer, int x, int y) {
        if (_free_idxs[layer] >= get_no_units_in_layer(layer)) {
            std::cerr << ">>> Lattice::set_pos_of_hidden_unit <<< no more free units in layer: " << layer << " of size: " << get_no_units_in_layer(layer) << std::endl;
            exit(EXIT_FAILURE);
        };

        int new_idx = _free_idxs[layer];
        _lookup_2[layer][x][y] = new_idx;
        _rlookup[layer][new_idx] = std::vector<int>({x,y});
        _free_idxs[layer] += 1;
    };
    void Lattice::set_pos_of_hidden_unit(int layer, int x, int y, int z) {
        if (_free_idxs[layer] >= get_no_units_in_layer(layer)) {
            std::cerr << ">>> Lattice::set_pos_of_hidden_unit <<< no more free units in layer: " << layer << " of size: " << get_no_units_in_layer(layer) << std::endl;
            exit(EXIT_FAILURE);
        };
        
        int new_idx = _free_idxs[layer];
        _lookup_3[layer][x][y][z] = new_idx;
        _rlookup[layer][new_idx] = std::vector<int>({x,y,z});
        _free_idxs[layer] += 1;
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
        _bias_dicts[layer][sp].push_back(bias);
        
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
        
        _o2_ixn_dicts[layer1][sp1][sp2].push_back(ixn);
        
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
        auto it = _bias_dicts.find(layer);
        double val = 0.0;
        if (it != _bias_dicts.end()) {
            auto it2 = it->second.find(sp);
            if (it2 != it->second.end()) {
                for (auto ixn: it2->second) {
                    val += ixn->get_val();
                };
            };
        };
        return val;
    };
    double Lattice::get_ixn_between_layers(int layer1, Sptr sp1, int layer2, Sptr sp2) const {
        if (layer2 != layer1+1 && layer2 != layer1-1) {
            std::cerr << ">>> Lattice::add_ixn_between_layers <<< layer2 != layer1 +- 1; instead layer_2 = " << layer2 << " and layer1 = " << layer1 << std::endl;
            exit(EXIT_FAILURE);
        };
        
        // Multiplier
        double mult = 1.0;
        auto itm = _o2_mults.find(layer1);
        if (itm != _o2_mults.end()) {
            auto itm2 = itm->second.find(layer2);
            if (itm2 != itm->second.end()) {
                mult = itm2->second;
            };
        };
        
        double val = 0.0;
        if (layer1 < layer2) {
            
            auto it = _o2_ixn_dicts.find(layer1);
            if (it != _o2_ixn_dicts.end()) {
                auto it2 = it->second.find(sp1);
                if (it2 != it->second.end()) {
                    auto it3 = it2->second.find(sp2);
                    if (it3 != it2->second.end()) {
                        for (auto ixn: it3->second) {
                            val += mult * ixn->get_val();
                        };
                    };
                };
            };
            
        } else {
            
            auto it = _o2_ixn_dicts.find(layer2);
            if (it != _o2_ixn_dicts.end()) {
                auto it2 = it->second.find(sp2);
                if (it2 != it->second.end()) {
                    auto it3 = it2->second.find(sp1);
                    if (it3 != it2->second.end()) {
                        for (auto ixn: it3->second) {
                            val += mult * ixn->get_val();
                        };
                    };
                };
            };
        };
        
        return val;
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
	void Lattice::set_empty_all_units() {
        for (auto layer=0; layer<_no_layers; layer++) {
            set_empty_all_units_in_layer(layer);
        };
	};
    void Lattice::set_empty_all_units_in_layer(int layer) {
        for (auto &sp_pair: _species_possible[layer]) {
            _latt[_i_markov_chain][layer][sp_pair.second].fill(0.0);
        };
    };

    // Random
    void Lattice::set_random_all_units(bool binary) {
        for (auto layer=0; layer<_no_layers; layer++) {
            set_random_all_units_in_layer(layer,binary);
        };
    };
    void Lattice::set_random_all_units_in_layer(int layer, bool binary) {
        std::cerr << ">>> Lattice::set_random_all_units_in_layer <<< not implemented!" << std::endl;
        exit(EXIT_FAILURE);

        if (binary) {
            // ...
        } else {
            // ...
        };
    };
    void Lattice::set_random_all_hidden_units(bool binary) {
        for (auto layer=1; layer<_no_layers; layer++) {
            set_random_all_units_in_layer(layer,binary);
        };
    };

    // Binarize
    void Lattice::binarize_all_units() {
        for (auto layer=0; layer<_no_layers; layer++) {
            binarize_all_units_in_layer(layer);
        };
    };
    void Lattice::binarize_all_units_in_layer(int layer) {
        // ...
        std::cerr << ">>> Lattice::binarize_all_units_in_layer <<< not implemented!" << std::endl;
        exit(EXIT_FAILURE);
    };

	/********************
	Write/read Latt to a file
	********************/

	void Lattice::write_layer_to_file(int layer, std::string fname, bool binary) const
	{
		std::ofstream f;
		f.open (fname);
		if (!f.is_open()) { // make sure we found it
			std::cerr << ">>> Lattice::write_layer_to_file <<< Error: could not open file: " << fname << " for writing" << std::endl;
			exit(EXIT_FAILURE);
		};
        
        // Go through all units
        int no_units = get_no_units_in_layer(layer);
        auto slatt = _latt.at(_i_markov_chain).at(layer);
        auto sps = _species_possible.at(layer);
        std::vector<int> pos;
        for (auto unit=0; unit<no_units; unit++) {
            // Write pos
            pos = _look_up_pos(layer,unit);
            for (auto const &x: pos) {
                f << x << " ";
            };
            
            // Go through species possible
            if (binary) {
                for (auto sp: sps) {
                    if (slatt.at(sp.second).at(unit) == 1.0) {
                        // Write species
                        f << sp.first << "\n";
                        break;
                    };
                };
            } else {
                for (auto sp: sps) {
                    // Write species
                    f << sp.first << " ";
                };
                f << "\n";
            };
        };
        
		f.close();
	};

	void Lattice::read_layer_from_file(int layer, std::string fname, bool binary)
	{
		// Clear the layer
        set_empty_all_units_in_layer(layer);

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
                for (auto i=0; i<_species_possible[layer].size(); i++) {
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
                for (auto i=0; i<_species_possible[layer].size(); i++) {
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
                for (auto i=0; i<_species_possible[layer].size(); i++) {
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
                auto it = _species_possible[layer].find(pr.first);
                if (it == _species_possible[layer].end()) {
                    std::cerr << ">>> Lattice::read_from_file <<< Error: could not find species (binary): " << pr.first << std::endl;
                    exit(EXIT_FAILURE);
                };
                for (auto const &site: pr.second) {
                    _latt[_i_markov_chain][layer][it->second][site] = 1.0;
                };
            };
        } else {
            // Prob mode
            for (auto const &pr: occs_to_write_prob) {
                // find the species
                auto it = _species_possible[layer].find(pr.first);
                if (it == _species_possible[layer].end()) {
                    std::cerr << ">>> Lattice::read_from_file <<< Error: could not find species (prob): " << pr.first << std::endl;
                    exit(EXIT_FAILURE);
                };
                for (auto const &site_pr: pr.second) {
                    _latt[_i_markov_chain][layer][it->second][site_pr.first] = site_pr.second;
                };
            };
        };
        
		f.close();
	};

	/********************
	Activate layer
	********************/

    // Activations
    void Lattice::_reset_activations(int layer) {
        for (auto &sp_pr: _species_possible.at(layer)) {
            _latt_act[_i_markov_chain][layer][sp_pr.second].fill(0.0);
        };
    };
    void Lattice::_calculate_activations(int layer, int given_layer) {
        // std::cout << "_calculate_activations: layer = " << layer << " given layer = " << given_layer << std::endl;
        int no_units = get_no_units_in_layer(layer);
        for (auto &sp_pr: _species_possible.at(layer)) {
            // std::cout << arma::size( get_bias_in_layer(layer, sp_pr.second) * arma::vec(no_units,arma::fill::ones) ) << std::endl;
            // std::cout << arma::size(_latt_act[_i_markov_chain][layer][sp_pr.second]) << std::endl;
            
            _latt_act[_i_markov_chain][layer][sp_pr.second] += get_bias_in_layer(layer, sp_pr.second) * arma::vec(no_units,arma::fill::ones);
            for (auto &given_sp_pr: _species_possible.at(given_layer)) {
                if (given_layer == layer-1) {
                    // Activate from below
                    _latt_act[_i_markov_chain][layer][sp_pr.second] += get_ixn_between_layers(given_layer, given_sp_pr.second, layer, sp_pr.second) * ( _adj[given_layer] * _latt_act[_i_markov_chain][given_layer][given_sp_pr.second] );
                } else if (given_layer == layer+1) {
                    // Activate from above
                    _latt_act[_i_markov_chain][layer][sp_pr.second] += get_ixn_between_layers(given_layer, given_sp_pr.second, layer, sp_pr.second) * ( _adj[layer].t() * _latt_act[_i_markov_chain][given_layer][given_sp_pr.second] );
                };
            };
        };
        // std::cout << "_calculate_activations: done!" << std::endl;
    };
    void Lattice::_convert_activations(int layer, bool binary) {
        
        // Convert activations to propensities via exp
        // Also calculate total propensity
        // Starts at 1.0 = exp(0) for empty
        int no_units = get_no_units_in_layer(layer);
        auto prop_tot = arma::vec(no_units,arma::fill::ones);
        for (auto &sp_pr: _species_possible.at(layer)) {
            _latt_act[_i_markov_chain][layer][sp_pr.second] = exp(_latt_act[_i_markov_chain][layer][sp_pr.second]);
            prop_tot += _latt_act[_i_markov_chain][layer][sp_pr.second];
        };
        
        // Divide by total to normalize
        // std::cout << "_convert_activations: layer: " << layer << " binary: " << binary << std::endl;
        for (auto &sp_pr: _species_possible.at(layer)) {
            _latt_act[_i_markov_chain][layer][sp_pr.second] /= prop_tot;
            // std::cout << sp_pr.first << " : " << _latt_act[_i_markov_chain][layer][sp_pr.second](0) << std::endl;
        };
        
        // Sample if binary
        if (binary) {
            // Random vec
            auto r = arma::vec(no_units,arma::fill::randu);
            
            // Convert to propensities
            auto running = arma::vec(no_units,arma::fill::zeros);
            for (auto &sp_pr: _species_possible.at(layer)) {
                running += _latt_act[_i_markov_chain][layer][sp_pr.second];
                _latt_act[_i_markov_chain][layer][sp_pr.second] = running;
            };
            
            // Evaluate
            // Initially: we are at zero, rand value is somewhere above
            // Next: we are at some propensity, rand value may be somewhere above or below
            // Finally: we are at max propensity = 1.0, rand value is below
            // flip vector:
            // 1 initially = rand is above current propensity
            // 0 = rand is below current propensity
            auto flip = arma::vec(no_units,arma::fill::ones);
            auto new_flip = flip;
            for (auto &sp_pr: _species_possible.at(layer)) {
                // New flip vector
                new_flip = 0.5*sign(r - _latt_act[_i_markov_chain][layer][sp_pr.second]) + 0.5;
                // flip -> new_flip
                // 1 -> 1 => not this species => (flip - new_flip) = 0
                // 1 -> 0 => yes this species => (flip - new_flip) = 1
                // 0 -> 0 => already assigned a previous species  => (flip - new_flip) = 0
                _latt_act[_i_markov_chain][layer][sp_pr.second] = flip - new_flip;
                
                // Old flip vector = new flip vector
                flip = new_flip;
            };
            // Finally: if still above the final propensity, it is the empty species! This means all are zero, so no further adjustment needed!
        };
    };
    
    // Prepare to activate a specific layer
    void Lattice::activate_layer_prepare(int layer, bool binary) {
        _reset_activations(layer);
        if (layer != 0) {
            _calculate_activations(layer, layer-1);
        };
        if (layer != _no_layers-1) {
            _calculate_activations(layer, layer+1);
        };
        _convert_activations(layer,binary);
    };
    void Lattice::activate_layer_prepare(int layer, int given_layer, bool binary) {
        if (given_layer != layer - 1 && given_layer != layer + 1) {
            std::cerr << ">>> Lattice::activate_layer_prepare <<< given layer must be +- layer, but instead layer = " << layer << " and given layer = " << given_layer << std::endl;
            exit(EXIT_FAILURE);
        };
        
        _reset_activations(layer);
        _calculate_activations(layer, given_layer);
        _convert_activations(layer,binary);
    };
    
    // Commit the activations
    void Lattice::activate_layer_committ(int layer) {
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
        
        _latt[_i_markov_chain][layer] = _latt_act[_i_markov_chain][layer];
    };
    
    // Activate a specific layer
    void Lattice::activate_layer(int layer, bool binary) {
        activate_layer_prepare(layer, binary);
        activate_layer_committ(layer);
    };
    void Lattice::activate_layer(int layer, int given_layer, bool binary) {
        activate_layer_prepare(layer, given_layer, binary);
        activate_layer_committ(layer);
    };
    
    // Variational inference
    void Lattice::variational_inference_hiddens() {
        
        // Prepare all layers
        for (auto layer=1; layer<_no_layers; layer++) {
            activate_layer_prepare(layer, false);
        };
        
        // Committ
        for (auto layer=1; layer<_no_layers; layer++) {
            activate_layer_committ(layer);
        };
    };
    
    // Sample
    void Lattice::sample(bool binary_visible, bool binary_hidden) {
        // Prepare all layers
        activate_layer_prepare(0, binary_visible);
        for (auto layer=1; layer<_no_layers; layer++) {
            activate_layer_prepare(layer, binary_hidden);
        };
        
        // Committ
        for (auto layer=0; layer<_no_layers; layer++) {
            activate_layer_committ(layer);
        };
    };
    
    // Make a pass activating upwards
    void Lattice::activate_upward_pass(bool binary_hidden) {
        for (auto layer=1; layer<_no_layers; layer++) {
            activate_layer(layer, layer-1, binary_hidden);
        };
    };
    
    /*
    void Lattice::sample_bm_up_v_to_h(bool binary_hidden, bool parallel) {
        for (auto &layer: _latt_h) {
            sample_layer(layer.first, layer.first-1, binary_hidden, parallel);
        };
    };

    // Variational inference in a BM
    void Lattice::sample_bm_variational_inference(bool parallel) {
        // Go through all hidden layers
        // Probabilistic mean field
        if (parallel) {
            
            // Parallel
            
            for (auto &it: _latt_h) {
                for (auto const &idx: _latt_h_idxs[it.first])
                {
                    it.second[idx]->prepare_sample(false); // false = probablistic
                };
            };
            for (auto &it: _latt_h) {
                for (auto const &idx: _latt_h_idxs[it.first])
                {
                    it.second[idx]->committ_sample();
                };
            };
        } else {
            
            // Not parallel
            
            for (auto &it: _latt_h) {
                
                // Shuffle
                std::random_shuffle ( _latt_h_idxs[it.first].begin(), _latt_h_idxs[it.first].end() );
                
                for (auto const &idx: _latt_h_idxs[it.first])
                {
                    it.second[idx]->prepare_sample(false); // false = probablistic
                    it.second[idx]->committ_sample();
                };
            };
        };
    };

    // Sample a specific layer
	void Lattice::sample_layer(int layer, bool binary, bool parallel) {

		if (layer == 0) {

			// Visible layer

			if (parallel) {

				// Parallel

				for (auto const &idx: _latt_v_idxs) 
				{
					_latt_v[idx]->prepare_sample(binary);
				};
				for (auto const &idx: _latt_v_idxs) 
				{
					_latt_v[idx]->committ_sample();
				};

			} else {

				// Not in parallel

				// Shuffle
				std::random_shuffle ( _latt_v_idxs.begin(), _latt_v_idxs.end() );		

				for (auto const &idx: _latt_v_idxs) 
				{
					_latt_v[idx]->prepare_sample(binary);
					_latt_v[idx]->committ_sample();
				};

			};

		} else {

			// Hidden layer
			auto it = _latt_h.find(layer);
			if (it == _latt_h.end()) {
				return;
			};

			if (parallel) {

				// Parallel

				for (auto const &idx: _latt_h_idxs[it->first]) 
				{
					it->second[idx]->prepare_sample(binary);
				};
				for (auto const &idx: _latt_h_idxs[it->first]) 
				{
					it->second[idx]->committ_sample();
				};

			} else {

				// Not in parallel

				std::random_shuffle ( _latt_h_idxs[it->first].begin(), _latt_h_idxs[it->first].end() );

				for (auto const &idx: _latt_h_idxs[it->first]) 
				{
					it->second[idx]->prepare_sample(binary);
					it->second[idx]->committ_sample();
				};
			};
		};
	};

	void Lattice::sample_layer(int layer, int given_layer, bool binary, bool parallel) {

		if (layer == 0) {

			// Visible layer

			if (given_layer != 1) {
				std::cerr << ">>> Lattice::sample_layer <<< Error: if sampling visible layer = 0, must specify given_layer = 1, instead of: " << given_layer << std::endl;
				exit(EXIT_FAILURE);
			};

			if (parallel) {

				// Parallel

				for (auto const &idx: _latt_v_idxs) 
				{
					_latt_v[idx]->prepare_sample(binary,given_layer);
				};
				for (auto const &idx: _latt_v_idxs) 
				{
					_latt_v[idx]->committ_sample();
				};

			} else {

				// Not in parallel

				// Shuffle
				std::random_shuffle ( _latt_v_idxs.begin(), _latt_v_idxs.end() );		

				for (auto const &idx: _latt_v_idxs) 
				{
					_latt_v[idx]->prepare_sample(binary,given_layer);
					_latt_v[idx]->committ_sample();
				};

			};

		} else {

			// Hidden layer

			if (given_layer != layer - 1 && given_layer != layer + 1) {
				std::cerr << ">>> Lattice::sample_layer <<< Error: sampling hidden layer: " << layer << " then given layer must be " << layer << " +- 1, instead of: " <<  given_layer << std::endl;
				exit(EXIT_FAILURE);
			};

			auto it = _latt_h.find(layer);
			if (it == _latt_h.end()) {
				return;
			};

			if (parallel) {

				// Parallel

				for (auto const &idx: _latt_h_idxs[it->first]) 
				{
					it->second[idx]->prepare_sample(binary,given_layer);
				};
				for (auto const &idx: _latt_h_idxs[it->first]) 
				{
					it->second[idx]->committ_sample();
				};

			} else {

				// Not in parallel

				std::random_shuffle ( _latt_h_idxs[it->first].begin(), _latt_h_idxs[it->first].end() );

				for (auto const &idx: _latt_h_idxs[it->first]) 
				{
					it->second[idx]->prepare_sample(binary,given_layer);
					it->second[idx]->committ_sample();
				};
			};
		};
	};

     */
     
    /********************
	Get counts
	********************/

	// 1 particle
	double Lattice::get_count_vis(Sptr &sp) const {
        return sum(_latt.at(_i_markov_chain).at(0).at(sp));
	};

    /*
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
    
    /********************
     Reap moments
     ********************/
    
    void Lattice::reap_moments(MomentType type, int i_sample) const {
        // Reset all ixns
        for (auto ixn: _all_ixns) {
            ixn->get_moment()->set_moment_sample(type, i_sample, 0.0);
        };
        
        // Reap biases
        double val;
        for (auto &bias_layer: _bias_dicts) {
            // Go through possible species in this layer
            for (auto &sp_pr: bias_layer.second) {
                // Get the moment
                val = arma::sum(_latt.at(_i_markov_chain).at(bias_layer.first).at(sp_pr.first));
        
                // Go through all ixns associated with this species
                for (auto ixn: sp_pr.second) {
                    // Set the moment
                    if (ixn->get_moment()) {
                        if (type == MomentType::ASLEEP || !ixn->get_moment()->get_is_awake_moment_fixed()) {
                            ixn->get_moment()->increment_moment_sample(type, i_sample, val);
                        };
                    };
                };
            };
        };
        
        // Reap ixns
        for (auto &o2_ixn_layer: _o2_ixn_dicts) {
            // Go through possible species in this layer
            for (auto &sp_pr_1: o2_ixn_layer.second) {
                for (auto &sp_pr_2: sp_pr_1.second) {
                    // Get the moment
                    val = dot(_latt.at(_i_markov_chain).at(o2_ixn_layer.first).at(sp_pr_1.first), _adj.at(o2_ixn_layer.first) * _latt.at(_i_markov_chain).at(o2_ixn_layer.first + 1).at(sp_pr_2.first));
                    
                    // Go through all ixns associated with this species
                    for (auto ixn: sp_pr_2.second) {
                        // Set the moment
                        if (ixn->get_moment()) {
                            if (type == MomentType::ASLEEP || !ixn->get_moment()->get_is_awake_moment_fixed()) {
                                ixn->get_moment()->increment_moment_sample(type, i_sample, val);
                            };
                        };
                    };
                };
            };
        };
    };
};
