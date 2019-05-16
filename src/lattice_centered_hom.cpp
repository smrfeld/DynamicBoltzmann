#include "../include/dblz_bits/lattice_centered_hom.hpp"

// Other headers
#include "../include/dblz_bits/general.hpp"
#include "../include/dblz_bits/species.hpp"
#include "../include/dblz_bits/moment_diff.hpp"
#include "../include/dblz_bits/ixn_param.hpp"
#include "../include/dblz_bits/fname.hpp"
#include "../include/dblz_bits/center.hpp"

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
	LatticeCenteredHom
	****************************************/

	/********************
	Constructor
	********************/

    LatticeCenteredHom::LatticeCenteredHom(int no_dims, int box_length, std::vector<Sptr> species_visible, std::vector<Cptr> centers) : Lattice(no_dims,box_length,species_visible,false) {
        
        // Visible layer
        add_layer(0, box_length, species_visible, centers, 0);
    };
    LatticeCenteredHom::LatticeCenteredHom(const LatticeCenteredHom& other) : Lattice(other) {
		_copy(other);
	};
    LatticeCenteredHom::LatticeCenteredHom(LatticeCenteredHom&& other) : Lattice(std::move(other)) {
		_move(other);
	};
	LatticeCenteredHom& LatticeCenteredHom::operator=(const LatticeCenteredHom& other) {
		if (this != &other) {
			_clean_up();
            Lattice::operator=(other);
			_copy(other);
		};
		return *this;
	};
	LatticeCenteredHom& LatticeCenteredHom::operator=(LatticeCenteredHom&& other) {
		if (this != &other) {
			_clean_up();
            Lattice::operator=(std::move(other));
			_move(other);
		};
		return *this;
	};
	LatticeCenteredHom::~LatticeCenteredHom() {
		_clean_up();
	};

	void LatticeCenteredHom::_clean_up() {};
	void LatticeCenteredHom::_copy(const LatticeCenteredHom& other) {
        _conn_mults = other._conn_mults;
        _centers = other._centers;
    };
	void LatticeCenteredHom::_move(LatticeCenteredHom& other) {
        _conn_mults = other._conn_mults;
        _centers = other._centers;
    
		// Reset other
        other._conn_mults.clear();
        other._centers.clear();
	};
    
    /********************
     Add a layer
     ********************/
    
    // This does nothing
    void LatticeCenteredHom::add_layer(int layer, int box_length, std::vector<Sptr> species) {
        std::cout << ">>> LatticeCenteredHom::add_layer <<< Error: this function does nothing without passing the centers." << std::endl;
        exit(EXIT_FAILURE);
    };

    // This is the correct function
    void LatticeCenteredHom::add_layer(int layer, int box_length, std::vector<Sptr> species, std::vector<Cptr> centers, int conn_multiplicity_to_layer_below) {
        
        if (layer != _no_layers) {
            std::cerr << ">>> LatticeCenteredHom::add_layer <<< error: next layer must be: " << _no_layers << " not: " << layer << std::endl;
            exit(EXIT_FAILURE);
        };
        
        // Conn mult
        if (layer != 0) {
            _conn_mults[layer-1][layer] = conn_multiplicity_to_layer_below;
            _conn_mults[layer][layer-1] = conn_multiplicity_to_layer_below;
        };

        // Increment no layers
        _no_layers += 1;
        
        // No units
        int no_units = pow(box_length,get_no_dims());
        _no_units_per_layer[layer] = no_units;

        // Persistent data structures
        _pst_prop[layer] = arma::vec(no_units);
        _pst_r[layer] = arma::vec(no_units);
        _pst_sign_of_r[layer] = arma::vec(no_units);
        _pst_sign_of_r_new[layer] = arma::vec(no_units);
        
        // Add zeros
        std::vector<MCType> chains({MCType::AWAKE,MCType::ASLEEP});
        for (auto const &chain: chains) {
            for (auto i_chain=0; i_chain<get_no_markov_chains(chain); i_chain++) {
                for (auto sp: species) {
                    _mc_chains[chain][i_chain][layer][sp] = arma::vec(no_units,arma::fill::zeros);
                    _mc_chains_act[chain][i_chain][layer][sp] = arma::vec(no_units,arma::fill::zeros);
                };
            };
        };
        
        // Lookups
        int ctr=0;
        if (get_no_dims() == 1) {
            for (auto x=1; x<=box_length; x++) {
                _lookup_1[layer][x] = ctr;
                _rlookup[layer][ctr] = std::vector<int>({x});
                ctr++;
            };
        } else if (get_no_dims() == 2) {
            for (auto x=1; x<=box_length; x++) {
                for (auto y=1; y<=box_length; y++) {
                    _lookup_2[layer][x][y] = ctr;
                    _rlookup[layer][ctr] = std::vector<int>({x,y});
                    ctr++;
                };
            };
        } else if (get_no_dims() == 3) {
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

        // Center
        for (auto center: centers) {
            _centers[center->get_layer()][center->get_species()] = center;
        };
    };
    
    // ***************
    // MARK: - Centers
    // ***************
    
    void LatticeCenteredHom::clear_all_centers() {
        _centers.clear();
    };
    void LatticeCenteredHom::set_centers_in_layer(int layer, std::vector<Cptr> centers) {
        _centers[layer].clear();
        for (auto center: centers) {
            _centers[layer][center->get_species()] = center;
        };
    };

    // ***************
    // MARK: - Reap moments
    // ***************
    
    void LatticeCenteredHom::reap_ixn_moment_diffs_and_slide_centers(double sliding_factor, bool calculate_offset) {
        
        // Reset all moments
        std::shared_ptr<MomentDiff> moment;
        for (auto ixn: _all_ixns) {
            moment = ixn->get_moment_diff();
            if (!moment->get_is_awake_moment_fixed()) {
                moment->reset_moment(MCType::AWAKE);
            };
            moment->reset_moment(MCType::ASLEEP);
        };
        
        // Reap biases
        // All biases
        int layer;
        Sptr sp;
        for (auto &bias_layer: _bias_dict) {
            layer = bias_layer.first;
            for (auto &sp_pr: bias_layer.second) {
                sp = sp_pr.first;
                moment = sp_pr.second->get_moment_diff();
                
                // Awake phase
                
                if (!moment->get_is_awake_moment_fixed()) {
                    for (auto i_chain=0; i_chain<get_no_markov_chains(MCType::AWAKE); i_chain++) {
                        moment->increment_moment(MCType::AWAKE, arma::accu(_mc_chains.at(MCType::AWAKE).at(i_chain).at(layer).at(sp)) / get_no_markov_chains(MCType::AWAKE));
                    };
                };
                
                // Asleep phase
                
                for (auto i_chain=0; i_chain<get_no_markov_chains(MCType::ASLEEP); i_chain++) {
                    moment->increment_moment(MCType::ASLEEP, arma::accu(_mc_chains.at(MCType::ASLEEP).at(i_chain).at(layer).at(sp)) / get_no_markov_chains(MCType::ASLEEP));
                };
            };
        };
        
        // Calculate the centers from the awake moment, and then slide
        int no_units;
        for (auto layer=0; layer<_no_layers; layer++) {
            no_units = get_no_units_in_layer(layer);
            
            // Determine means
            for (auto sp: _species_possible_vec.at(layer)) {
                // Set new val for center
                _centers.at(layer).at(sp)->set_val_new(_bias_dict.at(layer).at(sp)->get_moment_diff()->get_moment(MCType::AWAKE) / no_units);
                
                // Slide
                _centers.at(layer).at(sp)->slide(sliding_factor);
            };
        };
        
        // Reap ixns
        int layer1, layer2;
        Sptr sp1, sp2;
        arma::sp_mat::iterator mit, mit_end;
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
                        moment = sp_pr_2.second->get_moment_diff();
                        
                        // Awake
                        if (!moment->get_is_awake_moment_fixed()) {
                            for (auto i_chain=0; i_chain<get_no_markov_chains(MCType::AWAKE); i_chain++) {
                                mit = _adj.at(layer1).at(layer2).begin();
                                mit_end = _adj.at(layer1).at(layer2).end();
                                for(; mit != mit_end; ++mit) {
                                    moment->increment_moment(MCType::AWAKE, ( _mc_chains.at(MCType::AWAKE).at(i_chain).at(layer2).at(sp2)(mit.row()) - _centers.at(layer2).at(sp2)->get_val() ) * ( _mc_chains.at(MCType::AWAKE).at(i_chain).at(layer1).at(sp1)(mit.col()) - _centers.at(layer1).at(sp1)->get_val() ) / get_no_markov_chains(MCType::AWAKE) );
                                };
                            };
                        };
                        
                        // Asleep phase
                        for (auto i_chain=0; i_chain<get_no_markov_chains(MCType::ASLEEP); i_chain++) {
                            mit = _adj.at(layer1).at(layer2).begin();
                            mit_end = _adj.at(layer1).at(layer2).end();
                            for(; mit != mit_end; ++mit) {
                                moment->increment_moment(MCType::ASLEEP, ( _mc_chains.at(MCType::ASLEEP).at(i_chain).at(layer2).at(sp2)(mit.row()) - _centers.at(layer2).at(sp2)->get_val() ) * ( _mc_chains.at(MCType::ASLEEP).at(i_chain).at(layer1).at(sp1)(mit.col()) - _centers.at(layer1).at(sp1)->get_val() ) / get_no_markov_chains(MCType::ASLEEP) );
                            };
                        };
                    };
                };
            };
        };
        
        // Offset bias for centering
        if (calculate_offset) {
            
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
                            
                            offset += _conn_mults.at(layer-1).at(layer) * ixn->get_val() * _centers.at(layer-1).at(sp_below)->get_val();
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
                            offset += _conn_mults.at(layer).at(layer+1) * ixn->get_val() * _centers.at(layer+1).at(sp_above)->get_val();
                        };
                    };
                    
                    // Adjust bias in this layer
                    _bias_dict.at(layer).at(sp)->get_moment_diff()->set_moment_offset(offset);
                };
            };
        };
    };
    
    // ***************
    // MARK: - Get centers
    // ***************
    
    std::shared_ptr<Center> LatticeCenteredHom::get_center_for_species_in_layer(int layer, Sptr species) const {
        return _centers.at(layer).at(species);
    };
};
