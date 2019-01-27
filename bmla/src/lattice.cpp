#include "../include/bmla_bits/lattice.hpp"

// Other headers
#include "../include/bmla_bits/general.hpp"
#include "../include/bmla_bits/species.hpp"
#include "../include/bmla_bits/unit.hpp"
#include "../include/bmla_bits/connections.hpp"
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

        _IO_did_init = false;
        
		// Make a fully linked list of sites
        int i_chain = 0;
        int i_layer = 0;
		if (dim == 1) {
            for (auto sp: species_visible) {
                _latt[i_chain][i_layer][sp] = arma::vec(_box_length,0.0);
            };
        } else if (dim == 2) {
            for (auto sp: species_visible) {
                _latt[i_chain][i_layer][sp] = arma::vec(pow(_box_length,2),0.0);
            };
		} else if (dim == 3) {
            for (auto sp: species_visible) {
                _latt[i_chain][i_layer][sp] = arma::vec(pow(_box_length,3),0.0);
            };
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
        _latt = other._latt;
        _hlookup_1 = other._hlookup_1;
        _hlookup_2 = other._hlookup_2;
        _hlookup_3 = other._hlookup_3;
        _adj = other._adj;
        _species_possible = other._species_possible;
        _IO_did_init = other._IO_did_init;
    };
	void Lattice::_move(Lattice& other) {
        _no_dims = other._no_dims;
        _box_length = other._box_length;
        _no_markov_chains = other._no_markov_chains;
        _latt = other._latt;
        _hlookup_1 = other._hlookup_1;
        _hlookup_2 = other._hlookup_2;
        _hlookup_3 = other._hlookup_3;
        _adj = other._adj;
        _species_possible = other._species_possible;
        _IO_did_init = other._IO_did_init;

		// Reset other
		other._no_dims = 0;
		other._box_length = 0;
        other._no_markov_chains = 0;
        other._latt.clear();
		other._hlookup_1.clear();
		other._hlookup_2.clear();
		other._hlookup_3.clear();
        other._adj.clear();
        other._species_possible.clear();
        other._IO_did_init = false;
	};
    
    /****************************************
    PRIVATE METHODS
     ****************************************/

    // Lookup a site iterator from x,y,z
    int Lattice::_look_up_unit(int layer, int x) const {
        auto it = _hlookup_1.find(layer);
        if (it != _hlookup_1.end()) {
            auto it2 = it->second.find(x);
            if (it2 != it->second.end()) {
                return it2->second;
            };
        };
        
        std::cerr << ">>> Lattice::_look_up_unit <<< could not find layer: " << layer << " x: " << x << std::endl;
        exit(EXIT_FAILURE);
    };
    int Lattice::_look_up_unit(int layer, int x, int y) const {
        auto it = _hlookup_2.find(layer);
        if (it != _hlookup_2.end()) {
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
        auto it = _hlookup_3.find(layer);
        if (it != _hlookup_3.end()) {
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
    
    // Count helpers
    void Lattice::_get_count(double &count, Sptr &sp1, Sptr &sp2, const int &idx1, const int &idx2, const int &i_chain, bool binary, bool reversibly) const {
        
    };
    void Lattice::_get_count(double &count, Sptr &sp1, Sptr &sp2, Sptr &sp3, const int &idx1, const int &idx2, const int &idx3, const int &i_chain, bool binary, bool reversibly) const {
        
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
    
    /********************
     Add a layer
     ********************/
    
    void Lattice::add_layer(int layer, int no_units, std::vector<Sptr> species) {
        for (auto i_chain=0; i_chain<_no_markov_chains; i_chain++) {
            for (auto sp: species) {
                _latt[i_chain][layer][sp] = arma::vec(no_units,0.0);
            };
        };
        
        // Add species possible
        for (auto sp: species) {
            _species_possible[layer][sp->get_name()] = sp;
        };
    };
    
	/********************
	Helpers to setup all sites - Biases
	********************/

    // Biases
	void Lattice::set_bias_dict_all_units(std::shared_ptr<BiasDict> bias_dict) {
        for (auto layer=0; layer<_no_layers; layer++) {
            _bias_dicts[layer] = bias_dict;
        };
    };

    void Lattice::set_bias_dict_all_units_in_layer(int layer, std::shared_ptr<BiasDict> bias_dict) {
        _bias_dicts[layer] = bias_dict;
    };

    // Ixns
    void Lattice::set_ixn_dict_between_layers(int layer_1, int layer_2, std::shared_ptr<O2IxnDict> ixn_dict) {
        if (layer_2 != layer_1 + 1) {
            std::cerr << ">>> Lattice::set_ixn_dict_between_layers <<< only layer_2 = " << layer_2 << " = layer_1 + 1 = " << layer_1 + 1 << " is supported" << std::endl;
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
        _hlookup_1[layer][x] = new_idx;
    };
    void Lattice::add_hidden_unit(int layer, int x, int y) {
        int new_idx = _add_hidden_unit(layer);
        _hlookup_2[layer][x][y] = new_idx;
    };
    void Lattice::add_hidden_unit(int layer, int x, int y, int z) {
        int new_idx = _add_hidden_unit(layer);
        _hlookup_3[layer][x][y][z] = new_idx;
    };
    
	/********************
	Getters
	********************/

	int Lattice::get_no_units_v() { 
		return _latt_v.size(); 
	};
	int Lattice::get_no_units_h() { 
		int count=0;
		for (auto const &pr: _latt_h) {
			count += pr.second.size();
		};
		return count; 
	};

    /********************
     No Markov chains
     ********************/
    
    int Lattice::get_no_markov_chains() const {
        return _no_markov_chains;
    };
    void Lattice::set_no_markov_chains(int no_chains) {
        if (_no_markov_chains == no_chains) {
            return;
        };
        
        _no_markov_chains = no_chains;
        for (auto &s: _latt_v) {
            s->set_no_markov_chains(_no_markov_chains);
        };
        for (auto &pr: _latt_h) {
            for (auto &s: pr.second) {
                s->set_no_markov_chains(_no_markov_chains);
            };
        };
    };
    
    void Lattice::switch_to_markov_chain_no(int no) {
        for (auto &s: _latt_v) {
            s->switch_to_markov_chain_no(no);
        };
        for (auto &pr: _latt_h) {
            for (auto &s: pr.second) {
                s->switch_to_markov_chain_no(no);
            };
        };
    };
    void Lattice::switch_to_awake_statistics() {
        for (auto &s: _latt_v) {
            s->switch_to_awake_statistics();
        };
        for (auto &pr: _latt_h) {
            for (auto &s: pr.second) {
                s->switch_to_awake_statistics();
            };
        };
    };
    
	/********************
	Apply funcs to all units
	********************/

	// Clear the lattice
	void Lattice::all_units_v_set_empty() {
		for (auto &s: _latt_v) {
			s->set_empty();
		};
	};
	void Lattice::all_units_h_set_empty() {
		for (auto &pr: _latt_h) {
			for (auto &s: pr.second) {
				s->set_empty();
			};
		};
	};
	void Lattice::all_units_set_empty() {
		all_units_v_set_empty();
		all_units_h_set_empty();
	};
	void Lattice::all_units_in_layer_set_empty(int layer) {
		if (layer == 0) {
			all_units_v_set_empty();
		} else {
			auto it = _latt_h.find(layer);
			for (auto &s: it->second) {
				s->set_empty();
			};
		};
	};

	// Random
    void Lattice::all_units_random(bool binary) {
        all_units_v_random(binary);
        all_units_h_random(binary);
    };
    void Lattice::all_units_h_random(bool binary) {
        for (auto &it: _latt_h) {
            for (auto &ptr: it.second) {
                ptr->set_occ_random(binary);
            };
        };
    };
	void Lattice::all_units_v_random(bool binary) {
		for (auto &ptr: _latt_v) {
			ptr->set_occ_random(binary);
		};
	};

	void Lattice::all_units_in_layer_random(int layer, bool binary) {
		if (layer == 0) {
			for (auto &ptr: _latt_v) {
				ptr->set_occ_random(binary);
			};
		} else {
			auto it = _latt_h.find(layer);
			if (it != _latt_h.end()) {
				for (auto &ptr: it->second) {
					ptr->set_occ_random(binary);
				};
			};
		};
	};

	// Binarize
	void Lattice::all_units_v_binarize() {
		for (auto &ptr: _latt_v) {
			ptr->binarize();
		};
	};
	void Lattice::all_units_h_binarize() {
		for (auto &pr: _latt_h) {
			for (auto &s: pr.second) {
				s->binarize();
			};
		};
	};

	/********************
	Write/read Latt to a file
	********************/

	void Lattice::write_layer_to_file(int layer, std::string fname, bool binary) 
	{
		std::ofstream f;
		f.open (fname);
		if (!f.is_open()) { // make sure we found it
			std::cerr << ">>> Lattice::write_layer_to_file <<< Error: could not open file: " << fname << " for writing" << std::endl;
			exit(EXIT_FAILURE);
		};

		if (layer == 0) {

			// Visible layer

			for (auto const &l: _latt_v) {
				auto nonzero_map = l->get_nonzero_occs();
				// Only write nonzero
				if (nonzero_map.size() != 0) {

					if (_no_dims == 1) {
						f << l->x();
					} else if (_no_dims == 2) {
						f << l->x() << " " << l->y();
					} else if (_no_dims == 3) {
						f << l->x() << " " << l->y() << " " << l->z();
					};

					if (binary) {

						// Ensure binary
						if (nonzero_map.size() != 1) {
							std::cerr << ">>> Lattice::write_to_file <<< Error: indicated to write binary mode, but site is not binary!" << std::endl;
							exit(EXIT_FAILURE);
						};

						// binary
						for (auto const &pr: nonzero_map) {
							f << " " << pr.first->get_name();
						};
					} else {
						// multiple
						for (auto const &pr: nonzero_map) {
							f << " " << pr.first->get_name() << " " << pr.second;
						};
					};
					f << "\n";
				};
			};

		} else {

			// Hidden layer

			auto it = _latt_h.find(layer);
			for (auto const &l: it->second) {
				auto nonzero_map = l->get_nonzero_occs();
				// Only write nonzero
				if (nonzero_map.size() != 0) {

					if (_no_dims == 1) {
						f << l->x();
					} else if (_no_dims == 2) {
						f << l->x() << " " << l->y();
					} else if (_no_dims == 3) {
						f << l->x() << " " << l->y() << " " << l->z();
					};

					if (binary) {

						// Ensure binary
						if (nonzero_map.size() != 1) {
							std::cerr << ">>> Lattice::write_to_file <<< Error: indicated to write binary mode, but site is not binary!" << std::endl;
							exit(EXIT_FAILURE);
						};

						// binary
						for (auto const &pr: nonzero_map) {
							f << " " << pr.first->get_name();
						};
					} else {
						// multiple
						for (auto const &pr: nonzero_map) {
							f << " " << pr.first->get_name() << " " << pr.second;
						};
					};
					f << "\n";
				};
			};
		};

		f.close();
	};

	void Lattice::init_file_reader(std::map<int,std::vector<Sptr>> layers_species_possible) {
		for (auto pr: layers_species_possible) {
			for (auto sp: pr.second) {
				_species_possible[pr.first][sp->get_name()] = sp;
			};
		};
        
		_IO_did_init = true;
	};
	void Lattice::read_layer_from_file(int layer, std::string fname, bool binary)
	{        
		if (!_IO_did_init) {
			std::cerr << ">>> Lattice::read_layer_from_file <<< Error: run init_file_reader first!" << std::endl;
			exit(EXIT_FAILURE);
		};

		// Clear the layer
		all_units_in_layer_set_empty(layer);

		// Open
		std::ifstream f;
		f.open(fname);
		if (!f.is_open()) { // make sure we found it
			std::cerr << ">>> Error: Lattice::read_layer_from_file <<< could not find file: " << fname << std::endl;
			exit(EXIT_FAILURE);
		};

		if (layer == 0) {

			// Visible

			std::map<std::string, std::vector<UnitVisible*>> occs_to_write_binary;
			std::map<std::string, std::vector<std::pair<UnitVisible*,double>>> occs_to_write_prob;

			std::string x="",y="",z="";
			std::string sp="";
			std::string line;
			std::istringstream iss;
			UnitVisible* s;
			std::string prob="";
			double prob_val;

			if (_no_dims == 1 && binary) {
				while (getline(f,line)) {
					if (line == "") { continue; };
					iss = std::istringstream(line);				
					iss >> x >> sp;
		    		s = _look_up_unit_v(atoi(x.c_str()));
		    		occs_to_write_binary[sp].push_back(s);
		    		// Reset
			    	sp=""; x="";
			    };
			} else if (_no_dims == 1 && !binary) {
				while (getline(f,line)) {
					if (line == "") { continue; };
					iss = std::istringstream(line);				
					iss >> x;
		    		s = _look_up_unit_v(atoi(x.c_str()));
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
		    		s = _look_up_unit_v(atoi(x.c_str()),atoi(y.c_str()));
		    		occs_to_write_binary[sp].push_back(s);
		    		// Reset
			    	sp=""; x=""; y="";
			    };
			} else if (_no_dims == 2 && !binary) {
				while (getline(f,line)) {
					if (line == "") { continue; };
					iss = std::istringstream(line);				
					iss >> x >> y;
		    		s = _look_up_unit_v(atoi(x.c_str()),atoi(y.c_str()));
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
			    	s = _look_up_unit_v(atoi(x.c_str()),atoi(y.c_str()),atoi(z.c_str()));
		    		occs_to_write_binary[sp].push_back(s);
		    		// Reset
			    	sp=""; x=""; y=""; z="";
			    };
			} else if (_no_dims == 3 && !binary) {
				while (getline(f,line)) {
					if (line == "") { continue; };
					iss = std::istringstream(line);				
					iss >> x >> y >> z;
			    	s = _look_up_unit_v(atoi(x.c_str()),atoi(y.c_str()),atoi(z.c_str()));
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
						std::cerr << ">>> Lattice::read_from_file <<< Error: could not find species: " << pr.first << std::endl;
						exit(EXIT_FAILURE);
					};
					for (auto const &site: pr.second) {
						site->set_occ(it->second,1.0);
					};
				};
			} else {
				// Prob mode
				for (auto const &pr: occs_to_write_prob) {
					// find the species
					auto it = _species_possible[layer].find(pr.first);
					if (it == _species_possible[layer].end()) {
						std::cerr << ">>> Lattice::read_from_file <<< Error: could not find species: " << pr.first << std::endl;
						exit(EXIT_FAILURE);
					};
					for (auto const &site_pr: pr.second) {
						site_pr.first->set_occ(it->second,site_pr.second);
					};
				};
			};

		} else {

			// Hidden

			std::map<std::string, std::vector<UnitHidden*>> occs_to_write_binary;
			std::map<std::string, std::vector<std::pair<UnitHidden*,double>>> occs_to_write_prob;

			std::string x="",y="",z="";
			std::string sp="";
			std::string line;
			std::istringstream iss;
			UnitHidden* s;
			std::string prob="";
			double prob_val;

			if (_no_dims == 1 && binary) {
				while (getline(f,line)) {
					if (line == "") { continue; };
					iss = std::istringstream(line);				
					iss >> x >> sp;
		    		s = _look_up_unit_h(layer,atoi(x.c_str()));
		    		occs_to_write_binary[sp].push_back(s);
		    		// Reset
			    	sp=""; x="";
			    };
			} else if (_no_dims == 1 && !binary) {
				while (getline(f,line)) {
					if (line == "") { continue; };
					iss = std::istringstream(line);				
					iss >> x;
		    		s = _look_up_unit_h(layer,atoi(x.c_str()));
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
		    		s = _look_up_unit_h(layer,atoi(x.c_str()),atoi(y.c_str()));
		    		occs_to_write_binary[sp].push_back(s);
		    		// Reset
			    	sp=""; x=""; y="";
			    };
			} else if (_no_dims == 2 && !binary) {
				while (getline(f,line)) {
					if (line == "") { continue; };
					iss = std::istringstream(line);				
					iss >> x >> y;
		    		s = _look_up_unit_h(layer,atoi(x.c_str()),atoi(y.c_str()));
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
			    	s = _look_up_unit_h(layer,atoi(x.c_str()),atoi(y.c_str()),atoi(z.c_str()));
		    		occs_to_write_binary[sp].push_back(s);
		    		// Reset
			    	sp=""; x=""; y=""; z="";
			    };
			} else if (_no_dims == 3 && !binary) {
				while (getline(f,line)) {
					if (line == "") { continue; };
					iss = std::istringstream(line);				
					iss >> x >> y >> z;
			    	s = _look_up_unit_h(layer,atoi(x.c_str()),atoi(y.c_str()),atoi(z.c_str()));
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
						site->set_occ(it->second,1.0);
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
						site_pr.first->set_occ(it->second,site_pr.second);
					};
				};
			};
		};

		f.close();

		// std::cout << "Read: " << fname << std::endl;
	};

	/********************
	Sample
	********************/

    // Sample an RBM down/up
	void Lattice::sample_rbm_down_h_to_v(bool binary_visible, bool parallel) {
		sample_layer(0, 1, binary_visible, parallel);
	};
	void Lattice::sample_rbm_up_v_to_h(bool binary_hidden, bool parallel) {
        sample_layer(1, 0, binary_hidden, parallel);
    };

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

	/********************
	Get counts
	********************/

	// 1 particle
	double Lattice::get_count(Sptr &sp) const {
		double count = 0.0;
		for (auto const &s: _latt_v) {
			count += s->get_occ(sp);
		};

		return count;
	};

	// 2 particle
	double Lattice::get_count(Sptr &sp1, Sptr &sp2, bool reversibly) const {
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
	double Lattice::get_count(Sptr &sp1, Sptr &sp2, Sptr &sp3, bool reversibly) const {
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
	double Lattice::get_count(Sptr &sp1, Sptr &sp2, Sptr &sp3, Sptr &sp4, bool reversibly) const {
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

	/****************************************
	Lattice - PRIVATE
	****************************************/

	/********************
	Lookup a site iterator from x,y,z
	********************/

	UnitVisible* Lattice::_look_up_unit_v(int x) const {
		if (x > _box_length || x < 1) {
			return nullptr;
		};

		// Figure out index in list
		int n = x-1;

		return _latt_v[n];
	};
	UnitVisible* Lattice::_look_up_unit_v(int x, int y) const {
		if (x > _box_length || x < 1 || y > _box_length || y < 1) {
			return nullptr;
		};

		// Figure out index in list
		int n = (x-1)*_box_length + y-1;

		return _latt_v[n];
	};
	UnitVisible* Lattice::_look_up_unit_v(int x, int y, int z) const {
		if (x > _box_length || x < 1 || y > _box_length || y < 1 || z > _box_length || z < 1) {
			return nullptr;
		};

		// Figure out index in list
		int n = (x-1)*_box_length*_box_length + (y-1)*_box_length + z-1;

		return _latt_v[n];
	};

	UnitHidden* Lattice::_look_up_unit_h(int layer, int x) const {
		auto it1 = _hlookup_1.find(layer);
		if (it1 == _hlookup_1.end()) {
			return nullptr;
		};
		auto it2 = it1->second.find(x);
		if (it2 == it1->second.end()) {
			return nullptr;
		};

		return it2->second;
	};
	UnitHidden* Lattice::_look_up_unit_h(int layer, int x, int y) const {
		auto it1 = _hlookup_2.find(layer);
		if (it1 == _hlookup_2.end()) {
			return nullptr;
		};
		auto it2 = it1->second.find(x);
		if (it2 == it1->second.end()) {
			return nullptr;
		};
		auto it3 = it2->second.find(y);
		if (it3 == it2->second.end()) {
			return nullptr;
		};

		return it3->second;
	};
	UnitHidden* Lattice::_look_up_unit_h(int layer, int x, int y, int z) const {
		auto it1 = _hlookup_3.find(layer);
		if (it1 == _hlookup_3.end()) {
			return nullptr;
		};
		auto it2 = it1->second.find(x);
		if (it2 == it1->second.end()) {
			return nullptr;
		};
		auto it3 = it2->second.find(y);
		if (it3 == it2->second.end()) {
			return nullptr;
		};
		auto it4 = it3->second.find(z);
		if (it4 == it3->second.end()) {
			return nullptr;
		};

		return it4->second;
	};

	/********************
	Check dim
	********************/

	void Lattice::_check_dim(int dim) const {
		if (_no_dims != dim) {
			std::cerr << ">>> Error: Lattice::_check_dim <<< dim is: " << _no_dims << " but requested: " << dim << std::endl;
			exit(EXIT_FAILURE);
		};	
	};
};
