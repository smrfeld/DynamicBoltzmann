#include "../include/dblz_bits/lattice_alternating_binary.hpp"

// Other headers
#include "../include/dblz_bits/general.hpp"
#include "../include/dblz_bits/species.hpp"
#include "../include/dblz_bits/moment_diff.hpp"
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
	LatticeAlternatingBinary
	****************************************/

	/********************
	Constructor
	********************/

    LatticeAlternatingBinary::LatticeAlternatingBinary(int no_dims, int box_length, std::vector<Sptr> species_visible) : Lattice(no_dims, box_length, species_visible, true) {
    };
    LatticeAlternatingBinary::LatticeAlternatingBinary(const LatticeAlternatingBinary& other) : Lattice(other) {
		_copy(other);
	};
    LatticeAlternatingBinary::LatticeAlternatingBinary(LatticeAlternatingBinary&& other) : Lattice(std::move(other)) {
		_move(other);
	};
	LatticeAlternatingBinary& LatticeAlternatingBinary::operator=(const LatticeAlternatingBinary& other) {
		if (this != &other) {
			_clean_up();
            Lattice::operator=(other);
			_copy(other);
		};
		return *this;
	};
	LatticeAlternatingBinary& LatticeAlternatingBinary::operator=(LatticeAlternatingBinary&& other) {
		if (this != &other) {
			_clean_up();
            Lattice::operator=(std::move(other));
			_move(other);
		};
		return *this;
	};
	LatticeAlternatingBinary::~LatticeAlternatingBinary() {
		_clean_up();
	};

	void LatticeAlternatingBinary::_clean_up() {
        // Nothing...
	};
	void LatticeAlternatingBinary::_copy(const LatticeAlternatingBinary& other) {
    };
	void LatticeAlternatingBinary::_move(LatticeAlternatingBinary& other) {
	};
    
	/********************
	Apply funcs to all units
	********************/

    // Random
    void LatticeAlternatingBinary::set_random_all_units_in_layer(MCType chain, int i_chain, int layer, bool binary) {
        // Random probs
        for (auto sp: _species_possible_vec.at(layer)) {
            _mc_chains[chain][i_chain][layer][sp].fill(arma::fill::randu);
        };
        
        if (binary) {
            binarize_all_units_in_layer(chain, i_chain, layer);
        };
    };

    // Binarize
    void LatticeAlternatingBinary::_binarize_all_units_in_layer(MCType chain, int i_chain, int layer, bool act) {
        
        if (layer % 2 == 0) {
        
            // Multinomial layer
            
            // Random vec
            _pst_r.at(layer).fill(arma::fill::randu);

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
            
        } else {
            
            // Binary layers (seperate)
            
            // _latt or _latt_act?
            if (act) {
                
                for (auto sp: _species_possible_vec.at(layer)) {
                    
                    // Random vec [0,1]
                    _pst_r.at(layer).fill(arma::fill::randu);

                    // if:
                    // random r < act => occupied = 1
                    // r - act < 0 => sign(r-act) = -1 => 1
                    // else:
                    // random r > act => unoccupied = 0
                    // r - act > 0 => sign(r-act) = 1 => 0
                    // therefore:
                    // occupation = 0.5 * (1 - sign(r-act))
                    _mc_chains_act[chain][i_chain][layer][sp] = 0.5 * (1.0 - sign(_pst_r.at(layer) - _mc_chains_act.at(chain).at(i_chain).at(layer).at(sp)));
                };
                
            } else {
                
                for (auto sp: _species_possible_vec.at(layer)) {
                    _pst_r.at(layer).fill(arma::fill::randu);
                    _mc_chains[chain][i_chain][layer][sp] = 0.5 * (1.0 - sign(_pst_r.at(layer) - _mc_chains.at(chain).at(i_chain).at(layer).at(sp)));
                };

            };
        };
    };
    void LatticeAlternatingBinary::binarize_all_units_in_layer(MCType chain, int i_chain, int layer) {
        _binarize_all_units_in_layer(chain,i_chain,layer,false);
    };

	/********************
	Write/read Latt to a file
	********************/

	void LatticeAlternatingBinary::write_layer_to_file(MCType chain, int i_chain, int layer, std::string fname, bool binary) const
	{
        if (layer % 2 == 1) {
            std::cerr << ">>> LatticeAlternatingBinary::write_layer_to_file <<< Error: format unspecified here as of yet" << std::endl;
            exit(EXIT_FAILURE);
        };
        
		std::ofstream f;
		f.open (fname);
		if (!f.is_open()) { // make sure we found it
			std::cerr << ">>> LatticeAlternatingBinary::write_layer_to_file <<< Error: could not open file: " << fname << " for writing" << std::endl;
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

	void LatticeAlternatingBinary::read_layer_from_file(MCType chain, int i_chain, int layer, std::string fname, bool binary)
	{
        if (layer % 2 == 1) {
            std::cerr << ">>> LatticeAlternatingBinary::write_layer_to_file <<< Error: format unspecified here as of yet" << std::endl;
            exit(EXIT_FAILURE);
        };
        
		// Clear the layer
        set_empty_all_units_in_layer(chain, i_chain, layer);

		// Open
		std::ifstream f;
		f.open(fname);
		if (!f.is_open()) { // make sure we found it
			std::cerr << ">>> Error: LatticeAlternatingBinary::read_layer_from_file <<< could not find file: " << fname << std::endl;
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

        if (get_no_dims() == 1 && binary) {
            while (getline(f,line)) {
                if (line == "") { continue; };
                iss = std::istringstream(line);
                iss >> x >> sp;
                s = _look_up_unit(layer,atoi(x.c_str()));
                occs_to_write_binary[sp].push_back(s);
                // Reset
                sp=""; x="";
            };
        } else if (get_no_dims() == 1 && !binary) {
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
        } else if (get_no_dims() == 2 && binary) {
            while (getline(f,line)) {
                if (line == "") { continue; };
                iss = std::istringstream(line);
                iss >> x >> y >> sp;
                s = _look_up_unit(layer,atoi(x.c_str()),atoi(y.c_str()));
                occs_to_write_binary[sp].push_back(s);
                // Reset
                sp=""; x=""; y="";
            };
        } else if (get_no_dims() == 2 && !binary) {
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
        } else if (get_no_dims() == 3 && binary) {
            while (getline(f,line)) {
                if (line == "") { continue; };
                iss = std::istringstream(line);
                iss >> x >> y >> z >> sp;
                s = _look_up_unit(layer,atoi(x.c_str()),atoi(y.c_str()),atoi(z.c_str()));
                occs_to_write_binary[sp].push_back(s);
                // Reset
                sp=""; x=""; y=""; z="";
            };
        } else if (get_no_dims() == 3 && !binary) {
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
                    std::cerr << ">>> LatticeAlternatingBinary::read_from_file <<< Error: could not find species (binary): " << pr.first << std::endl;
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
                    std::cerr << ">>> LatticeAlternatingBinary::read_from_file <<< Error: could not find species (prob): " << pr.first << std::endl;
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
    

    // (2) Convert activations to probs
    void LatticeAlternatingBinary::activate_layer_convert_to_probs(MCType chain, int layer, bool binary) {
        
        if (layer % 2 == 0) {
            
            // Multinomial layer
            
            // All chains
            for (auto i_chain=0; i_chain<get_no_markov_chains(chain); i_chain++) {
            
                // Convert activations to propensities via exp
                // Also calculate total propensity
                _pst_prop.at(layer).fill(arma::fill::ones); // Starts at 1.0 = exp(0) for empty
                for (auto sp: _species_possible_vec.at(layer)) {
                    _mc_chains_act[chain][i_chain][layer][sp].transform( [](double val) { return (exp(val)); } );
                    _pst_prop[layer] += _mc_chains_act.at(chain).at(i_chain).at(layer).at(sp);
                };
                
                // Divide by total to normalize
                for (auto sp: _species_possible_vec.at(layer)) {
                    _mc_chains_act[chain][i_chain][layer][sp] /= _pst_prop.at(layer);
                };
            };
            
        } else {
            
            // Binary layer
            
            // Fill 1 for empty
            _pst_prop.at(layer).fill(arma::fill::ones);
            
            // All chains
            for (auto i_chain=0; i_chain<get_no_markov_chains(chain); i_chain++) {
                
                // Convert activations to propensities via exp
                for (auto sp: _species_possible_vec.at(layer)) {
                    _mc_chains_act[chain][i_chain][layer][sp].transform( [](double val) { return (exp(val)); } );
                };
                
                // Divide by total to normalize
                for (auto sp: _species_possible_vec.at(layer)) {
                    _mc_chains_act[chain][i_chain][layer][sp] /= _pst_prop.at(layer) + _mc_chains_act.at(chain).at(i_chain).at(layer).at(sp);
                };
            };
        };
            
        // Sample if binary
        if (binary) {
            for (auto i_chain=0; i_chain<get_no_markov_chains(chain); i_chain++) {
                _binarize_all_units_in_layer(chain,i_chain,layer,true);
            };
        };
    };
};
