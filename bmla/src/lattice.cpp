#include "../include/bmla_bits/lattice.hpp"

// Other headers
#include "../include/bmla_bits/general.hpp"
#include "../include/bmla_bits/species.hpp"
#include "../include/bmla_bits/ixn_dicts.hpp"
#include "../include/bmla_bits/moment.hpp"

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
        _no_layers = 1;
        _i_markov_chain = 0;
        
		// Visible layer
        int i_chain = 0;
        int i_layer = 0;
		if (dim == 1) {
            for (auto sp: species_visible) {
                _latt[i_chain][i_layer][sp] = arma::vec(_box_length,0.0);
            };
            for (auto x=1; x<=_box_length; x++) {
                _lookup_1[0][x] = x-1;
                _rlookup[0][x-1] = std::vector<int>({x});
            };
        } else if (dim == 2) {
            for (auto sp: species_visible) {
                _latt[i_chain][i_layer][sp] = arma::vec(pow(_box_length,2),0.0);
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
                _latt[i_chain][i_layer][sp] = arma::vec(pow(_box_length,3),0.0);
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
        _no_layers = other._no_layers;
        _i_markov_chain = other._i_markov_chain;
        
        _latt = other._latt;
        _lookup_1 = other._lookup_1;
        _lookup_2 = other._lookup_2;
        _lookup_3 = other._lookup_3;
        _rlookup = other._rlookup;
        _adj = other._adj;
        
        _bias_dicts = other._bias_dicts;
        _ixn_dicts = other._ixn_dicts;
        
        _species_possible = other._species_possible;
    };
	void Lattice::_move(Lattice& other) {
        _no_dims = other._no_dims;
        _box_length = other._box_length;
        _no_markov_chains = other._no_markov_chains;
        _no_layers = other._no_layers;
        _i_markov_chain = other._i_markov_chain;
        
        _latt = other._latt;
        _lookup_1 = other._lookup_1;
        _lookup_2 = other._lookup_2;
        _lookup_3 = other._lookup_3;
        _rlookup = other._rlookup;
        _adj = other._adj;
        
        _bias_dicts = other._bias_dicts;
        _ixn_dicts = other._ixn_dicts;
        
        _species_possible = other._species_possible;

		// Reset other
		other._no_dims = 0;
		other._box_length = 0;
        other._no_markov_chains = 0;
        other._no_layers = 0;
        other._i_markov_chain = 0;
        
        other._latt.clear();
		other._lookup_1.clear();
		other._lookup_2.clear();
		other._lookup_3.clear();
        other._rlookup.clear();
        other._adj.clear();
        
        other._bias_dicts.clear();
        other._ixn_dicts.clear();
        
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
    
    int Lattice::get_no_markov_chains() const {
        return _no_markov_chains;
    };
    void Lattice::set_no_markov_chains(int no_markov_chains) {
        _no_markov_chains = no_markov_chains;
        
        if (_latt.size() < _no_markov_chains) {
            for (auto i_chain=_latt.size(); i_chain < _no_markov_chains; i_chain++) {
                _latt[i_chain] = _latt[i_chain-1];
            };
        };
        if (_latt.size() > _no_markov_chains) {
            for (auto i_chain=_latt.size(); i_chain > _no_markov_chains; i_chain--) {
                _latt.erase(i_chain);
            };
        };
    };
    void Lattice::set_current_markov_chain(int i_markov_chain) {
        _i_markov_chain = i_markov_chain;
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
                _latt[i_chain][layer][sp] = arma::vec(no_units,0.0);
            };
        };
        
        // Add species possible
        for (auto sp: species) {
            _species_possible[layer][sp->get_name()] = sp;
        };
        
        // Add adjacency matrix
        int size_below = get_no_units_in_layer(layer-1);
        _adj[layer-1] = arma::mat(size_below,no_units,arma::fill::zeros);
    };
    
	/********************
	Helpers to setup all sites - Biases
	********************/

    // Biases
	void Lattice::set_bias_dict_all_units(std::shared_ptr<BiasDict> bias_dict) {
        for (auto layer=0; layer<_no_layers; layer++) {
            set_bias_dict_all_units_in_layer(layer, bias_dict);
        };
    };

    void Lattice::set_bias_dict_all_units_in_layer(int layer, std::shared_ptr<BiasDict> bias_dict) {
        _bias_dicts[layer] = bias_dict;
    };

    // Ixns
    void Lattice::set_ixn_dict_between_layers(int layer_1, int layer_2, std::shared_ptr<O2IxnDict> ixn_dict) {
        if (layer_2 != layer_1 + 1) {
            std::cerr << ">>> Lattice::set_ixn_dict_between_layers <<< only layer_2 = " << layer_2 << " = layer_1 + 1 = " << layer_1 + 1 << " is supported" << std::endl;
            exit(EXIT_FAILURE);
        };
        _ixn_dicts[layer_1] = ixn_dict;
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
     Add hidden units
     ********************/
    
    // Add hidden unit to layer
    int Lattice::_add_hidden_unit(int layer) {
        // Add unit (resize)
        int size_current = 0;
        for (auto i_chain=0; i_chain<_no_markov_chains; i_chain++) {
            for (auto &sp_pair: _species_possible[layer]) {
                size_current = _latt[i_chain][layer][sp_pair.second].size();
                _latt[i_chain][layer][sp_pair.second].resize(size_current + 1);
            };
        };
        
        // Resize adjacency matrix
        // The one below
        int size1=0, size2=0;
        if (layer != 0) {
            size1 = _adj[layer-1].n_rows;
            size2 = _adj[layer-1].n_cols;
            _adj[layer-1].resize(size1, size2+1);
        };
        // The one above
        if (layer != _no_layers-1) {
            size1 = _adj[layer].n_rows;
            size2 = _adj[layer].n_cols;
            _adj[layer].resize(size1+1, size2);
        };
        
        return size_current;
    };
    
    void Lattice::add_hidden_unit(int layer, int x) {
        int new_idx = _add_hidden_unit(layer);
        _lookup_1[layer][x] = new_idx;
        _rlookup[layer][new_idx] = std::vector<int>({x});
    };
    void Lattice::add_hidden_unit(int layer, int x, int y) {
        int new_idx = _add_hidden_unit(layer);
        _lookup_2[layer][x][y] = new_idx;
        _rlookup[layer][new_idx] = std::vector<int>({x,y});
    };
    void Lattice::add_hidden_unit(int layer, int x, int y, int z) {
        int new_idx = _add_hidden_unit(layer);
        _lookup_3[layer][x][y][z] = new_idx;
        _rlookup[layer][new_idx] = std::vector<int>({x,y,z});
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
        if (binary) {
            // ...
        } else {
            // ...
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

    // Activate a specific layer
    void Lattice::activate_layer(int layer, bool binary) {
        
    };
    void Lattice::activate_layer(int layer, int given_layer, bool binary) {
        
    };
    
    /*
    // Sample BM (NOT layer-wise)
    void Lattice::sample_bm(bool binary_visible, bool binary_hidden, bool parallel) {
        
        if (parallel) {
            
            // Parallel
            
            // Visible
            for (auto const &idx: _latt_v_idxs)
            {
                _latt_v[idx]->prepare_sample(binary_visible);
            };
            
            // Hiddens
            for (auto &it: _latt_h) {
                for (auto const &idx: _latt_h_idxs[it.first]) {
                    it.second[idx]->prepare_sample(binary_hidden);
                };
            };
            
            // Committ
            for (auto const &idx: _latt_v_idxs)
            {
                _latt_v[idx]->committ_sample();
            };
            for (auto &it: _latt_h) {
                for (auto const &idx: _latt_h_idxs[it.first]) {
                    it.second[idx]->committ_sample();
                };
            };

        } else {
            
            // Not parallel
            
            // Shuffle visible
            std::random_shuffle ( _latt_v_idxs.begin(), _latt_v_idxs.end() );
            
            // Visible
            for (auto const &idx: _latt_v_idxs)
            {
                _latt_v[idx]->prepare_sample(binary_visible);
                _latt_v[idx]->committ_sample();
            };
            
            // Hiddens
            for (auto &it: _latt_h) {
                // Shuffle hidden
                std::random_shuffle ( _latt_h_idxs[it.first].begin(), _latt_h_idxs[it.first].end() );
                
                for (auto const &idx: _latt_h_idxs[it.first]) {
                    it.second[idx]->prepare_sample(binary_hidden);
                    it.second[idx]->committ_sample();
                };
            };
        };
    };
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
};
