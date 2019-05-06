#include "../include/dblz_bits/lattice_1d_fully_visible.hpp"

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
	Lattice1DFullyVisible
	****************************************/

	/********************
	Constructor
	********************/

    Lattice1DFullyVisible::Lattice1DFullyVisible(int box_length, std::vector<Sptr> species_visible, bool conn_2, bool conn_3) {
        
        _box_length = box_length;
        _conn_2 = conn_2;
        _conn_3 = conn_3;
        
        // Set no chains
        set_no_markov_chains(MCType::AWAKE, 1);
        set_no_markov_chains(MCType::ASLEEP, 1);
        
        // Visible layer
        
        // Add random
        for (auto &chain: _no_markov_chains) {
            for (auto i_chain=0; i_chain<chain.second; i_chain++) {
                for (auto sp: species_visible) {
                    _mc_chains[chain.first][i_chain][sp] = arma::vec(_box_length,arma::fill::randu);
                    _mc_chains_act[chain.first][i_chain][sp] = arma::vec(_box_length,arma::fill::randu);
                };
            };
        };
        
        // Lookups
        int ctr=0;
        for (auto x=1; x<=box_length; x++) {
            _lookup_1[x] = ctr;
            _rlookup[ctr] = x;
            ctr++;
        };
        
        // Add species possible
        for (auto sp: species_visible) {
            _species_possible_map[sp->get_name()] = sp;
            _species_possible_vec.push_back(sp);
        };
        
        // Persistent data structures
        _pst_prop = arma::vec(box_length);
        _pst_r = arma::vec(box_length);
        _pst_sign_of_r = arma::vec(box_length);
        _pst_sign_of_r_new = arma::vec(box_length);
    };
	Lattice1DFullyVisible::Lattice1DFullyVisible(const Lattice1DFullyVisible& other) {
		_copy(other);
	};
	Lattice1DFullyVisible::Lattice1DFullyVisible(Lattice1DFullyVisible&& other) {
		_move(other);
	};
	Lattice1DFullyVisible& Lattice1DFullyVisible::operator=(const Lattice1DFullyVisible& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	Lattice1DFullyVisible& Lattice1DFullyVisible::operator=(Lattice1DFullyVisible&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
	};
	Lattice1DFullyVisible::~Lattice1DFullyVisible() {
		_clean_up();
	};

	void Lattice1DFullyVisible::_clean_up() {
        // Nothing...
	};
	void Lattice1DFullyVisible::_copy(const Lattice1DFullyVisible& other) {
		_box_length = other._box_length;
        _no_markov_chains = other._no_markov_chains;
        
        _mc_chains = other._mc_chains;
        _mc_chains_act = other._mc_chains_act;
        
        _lookup_1 = other._lookup_1;
        _rlookup = other._rlookup;
        
        _species_possible_map = other._species_possible_map;
        _species_possible_vec = other._species_possible_vec;
        
        _conn_2 = other._conn_2;
        _conn_3 = other._conn_3;
        
        _all_ixns = other._all_ixns;
        _bias_dict = other._bias_dict;
        _o2_ixn_dict = other._o2_ixn_dict;
        _o3_ixn_dict = other._o3_ixn_dict;
        
        _pst_prop = other._pst_prop;
        _pst_r = other._pst_r;
        _pst_sign_of_r = other._pst_sign_of_r;
        _pst_sign_of_r_new = other._pst_sign_of_r_new;
    };
	void Lattice1DFullyVisible::_move(Lattice1DFullyVisible& other) {
        _box_length = other._box_length;
        _no_markov_chains = other._no_markov_chains;

        _mc_chains = other._mc_chains;
        _mc_chains_act = other._mc_chains_act;
        
        _lookup_1 = other._lookup_1;
        _rlookup = other._rlookup;
        
        _species_possible_map = other._species_possible_map;
        _species_possible_vec = other._species_possible_vec;
        
        _conn_2 = other._conn_2;
        _conn_3 = other._conn_3;

        _all_ixns = other._all_ixns;
        _bias_dict = other._bias_dict;
        _o2_ixn_dict = other._o2_ixn_dict;
        _o3_ixn_dict = other._o3_ixn_dict;

        _pst_prop = other._pst_prop;
        _pst_r = other._pst_r;
        _pst_sign_of_r = other._pst_sign_of_r;
        _pst_sign_of_r_new = other._pst_sign_of_r_new;
        
		// Reset other
		other._box_length = 0;
        other._no_markov_chains.clear();
        
        other._mc_chains.clear();
        other._mc_chains_act.clear();
        
		other._lookup_1.clear();
        other._rlookup.clear();
        
        other._species_possible_map.clear();
        other._species_possible_vec.clear();
        
        other._conn_2 = false;
        other._conn_3 = false;
        
        other._all_ixns.clear();
        other._bias_dict.clear();
        other._o2_ixn_dict.clear();
        other._o3_ixn_dict.clear();
        
        other._pst_prop.clear();
        other._pst_r.clear();
        other._pst_sign_of_r.clear();
        other._pst_sign_of_r_new.clear();
	};
    
    /****************************************
    PRIVATE METHODS
     ****************************************/

    // Lookup a site iterator from x,y,z
    int Lattice1DFullyVisible::_look_up_unit(int x) const {
        auto it = _lookup_1.find(x);
        if (it != _lookup_1.end()) {
            return it->second;
        };
        
        std::cerr << ">>> Lattice1DFullyVisible::_look_up_unit <<< could not find x: " << x << std::endl;
        exit(EXIT_FAILURE);
    };
    
    int Lattice1DFullyVisible::_look_up_pos(int idx) const {
        auto it = _rlookup.find(idx);
        if (it != _rlookup.end()) {
            return it->second;
        };

        return 0;
    };
    
    /****************************************
    PUBLIC METHODS
     ****************************************/

    /********************
     Getters
     ********************/
    
    int Lattice1DFullyVisible::get_box_length() const {
        return _box_length;
    };
    
    int Lattice1DFullyVisible::get_no_units() const {
        return _box_length;
    };
    
    /********************
     Markov chains
     ********************/
    
    int Lattice1DFullyVisible::get_no_markov_chains(MCType type) const {
        return _no_markov_chains.at(type);
    };

    void Lattice1DFullyVisible::set_no_markov_chains(MCType type, int no_markov_chains) {
        _no_markov_chains[type] = no_markov_chains;
        
        // Add chains
        if (_mc_chains[type].size() < _no_markov_chains[type]) {
            for (auto i_chain=_mc_chains[type].size(); i_chain < _no_markov_chains[type]; i_chain++) {
                for (auto sp: _species_possible_vec) {
                    _mc_chains[type][i_chain][sp] = arma::vec(_box_length,arma::fill::randu);
                    _mc_chains_act[type][i_chain][sp] = arma::vec(_box_length,arma::fill::randu);
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
	Helpers to setup all sites - Biases
	********************/

    // Biases
    void Lattice1DFullyVisible::set_bias(Sptr sp, Iptr bias) {
        _bias_dict[sp] = bias;
        
        _add_to_all_ixns_vec(bias);
    };

    // Ixns
    void Lattice1DFullyVisible::set_ixn_2(Sptr sp1, Sptr sp2, Iptr ixn) {
        // Add both ways
        _o2_ixn_dict[sp1][sp2] = ixn;
        _o2_ixn_dict[sp2][sp1] = ixn;

        // Add to all
        _add_to_all_ixns_vec(ixn);
    };
    void Lattice1DFullyVisible::set_ixn_3(Sptr sp1, Sptr sp2, Sptr sp3, Iptr ixn) {
        // Add both ways
        _o3_ixn_dict[sp1][sp2][sp3] = ixn;
        _o3_ixn_dict[sp3][sp2][sp1] = ixn;
        
        // Add to all
        _add_to_all_ixns_vec(ixn);
    };
    
    // Get ixns
    double Lattice1DFullyVisible::get_bias(Sptr sp) const {
        auto it = _bias_dict.find(sp);
        if (it != _bias_dict.end()) {
            return it->second->get_val();
        } else {
            return 0.0;
        };
    };
    double Lattice1DFullyVisible::get_ixn_2(Sptr sp1, Sptr sp2) const {
        auto it1 = _o2_ixn_dict.find(sp1);
        if (it1 != _o2_ixn_dict.end()) {
            auto it2 = it1->second.find(sp2);
            if (it2 != it1->second.end()) {
                return it2->second->get_val();
            };
        };
        return 0.0;
    };
    double Lattice1DFullyVisible::get_ixn_3(Sptr sp1, Sptr sp2, Sptr sp3) const {
        auto it1 = _o3_ixn_dict.find(sp1);
        if (it1 != _o3_ixn_dict.end()) {
            auto it2 = it1->second.find(sp2);
            if (it2 != it1->second.end()) {
                auto it3 = it2->second.find(sp3);
                if (it3 != it2->second.end()) {
                    return it3->second->get_val();
                };
            };
        };
        return 0.0;
    };

    // Get all ixns
    const std::vector<Iptr>& Lattice1DFullyVisible::get_all_ixn_params() const {
        return _all_ixns;
    };
    
	/********************
	Apply funcs to all units
	********************/

	// Clear the Lattice1DFullyVisible
	void Lattice1DFullyVisible::set_empty_all_units(MCType chain, int i_chain) {
        for (auto sp: _species_possible_vec) {
            _mc_chains[chain][i_chain][sp].fill(0.0);
        };
	};

    // Random
    void Lattice1DFullyVisible::set_random_all_units(MCType chain, int i_chain, bool binary) {
        // Random probs
        for (auto sp: _species_possible_vec) {
            _mc_chains[chain][i_chain][sp].fill(arma::fill::randu);
        };
        
        if (binary) {
            binarize_all_units(chain, i_chain, true);
        };
    };

    // Binarize
    int _sign(double val) {
        return (0.0 < val) - (val < 0.0);
    }
    void Lattice1DFullyVisible::binarize_unit(MCType chain, int i_chain, int x, bool act) {
        // Random vec
        double _pst_r = randD(0.0,1.0);
        
        // Helpers
        
        // sign of r vector:
        // 1 initially
        // 1 -> 1 => not this species
        // 1 -> -1 => this species
        // -1 -> -1 => not this species
        double _pst_sign_of_r = 1;
        double _pst_sign_of_r_new = 1;
        
        // _latt or _latt_act?
        if (act) {
            
            // Evaluate
            // Initially: we (flip) are at zero, rand value is somewhere above
            // Next: we (flip) are at some propensity, rand value may be somewhere above or below
            // Finally: we (flip) are at max propensity = 1.0, rand value is below
            for (auto sp: _species_possible_vec) {
                // Subtract prob from r
                _pst_r -= _mc_chains_act[chain][i_chain][sp][x];
                
                // New flip vector
                _pst_sign_of_r_new = _sign(_pst_r);
                
                // Difference
                // sign_of_r -> sign_of_r_new
                // 1 -> 1 => not this species [0]
                // 1 -> -1 => this species [1]
                // -1 -> -1 => not this species [0]
                // 0.5 * (sign_of_r - sign_of_r_new) encodes the [brackets]
                
                _mc_chains_act[chain][i_chain][sp][x] = 0.5 * (_pst_sign_of_r - _pst_sign_of_r_new);
                
                // Old flip vector = new flip vector
                _pst_sign_of_r = _pst_sign_of_r_new;
            };
            // Finally: if still above the final propensity, it is the empty species! This means all are zero, so no further adjustment needed!
            
        } else {
            
            for (auto sp: _species_possible_vec) {
                _pst_r -= _mc_chains[chain][i_chain][sp][x];
                _pst_sign_of_r_new = _sign(_pst_r);
                _mc_chains[chain][i_chain][sp][x] = 0.5 * (_pst_sign_of_r - _pst_sign_of_r_new);
                _pst_sign_of_r = _pst_sign_of_r_new;
            };
        };
    };
    void Lattice1DFullyVisible::binarize_all_units(MCType chain, int i_chain, bool act) {
        // Random vec
        _pst_r.fill(arma::fill::randu);
        
        // Helpers
        
        // sign of r vector:
        // 1 initially
        // 1 -> 1 => not this species
        // 1 -> -1 => this species
        // -1 -> -1 => not this species
        _pst_sign_of_r.fill(arma::fill::ones);
        _pst_sign_of_r_new.fill(arma::fill::ones);
        
        // _latt or _latt_act?
        if (act) {
            
            // Evaluate
            // Initially: we (flip) are at zero, rand value is somewhere above
            // Next: we (flip) are at some propensity, rand value may be somewhere above or below
            // Finally: we (flip) are at max propensity = 1.0, rand value is below
            for (auto sp: _species_possible_vec) {
                // Subtract prob from r
                _pst_r -= _mc_chains_act[chain][i_chain][sp];
                
                // New flip vector
                _pst_sign_of_r_new = sign(_pst_r);
                
                // Difference
                // sign_of_r -> sign_of_r_new
                // 1 -> 1 => not this species [0]
                // 1 -> -1 => this species [1]
                // -1 -> -1 => not this species [0]
                // 0.5 * (sign_of_r - sign_of_r_new) encodes the [brackets]
                
                _mc_chains_act[chain][i_chain][sp] = 0.5 * (_pst_sign_of_r - _pst_sign_of_r_new);
                
                // Old flip vector = new flip vector
                _pst_sign_of_r = _pst_sign_of_r_new;
            };
            // Finally: if still above the final propensity, it is the empty species! This means all are zero, so no further adjustment needed!
            
        } else {
            
            for (auto sp: _species_possible_vec) {
                _pst_r -= _mc_chains[chain][i_chain][sp];
                _pst_sign_of_r_new = sign(_pst_r);
                _mc_chains[chain][i_chain][sp] = 0.5 * (_pst_sign_of_r - _pst_sign_of_r_new);
                _pst_sign_of_r = _pst_sign_of_r_new;
            };
        };
    };

	/********************
	Write/read Latt to a file
	********************/

	void Lattice1DFullyVisible::write_to_file(MCType chain, int i_chain, std::string fname, bool binary) const
	{
		std::ofstream f;
		f.open (fname);
		if (!f.is_open()) { // make sure we found it
			std::cerr << ">>> Lattice1DFullyVisible::write_layer_to_file <<< Error: could not open file: " << fname << " for writing" << std::endl;
			exit(EXIT_FAILURE);
		};
        
        // Go through all units
        auto slatt = _mc_chains.at(chain).at(i_chain);
        auto sps = _species_possible_vec;
        if (binary) {

            // Binary

            // Go through all units
            for (auto unit=0; unit<_box_length; unit++) {
                
                // Go through species possible
                for (auto sp: sps) {
                    if (abs(slatt.at(sp)(unit) - 1.0) < 1.0e-5) {
                        // Write pos
                        f << _look_up_pos(unit) << " ";
                        
                        // Write species
                        f << sp->get_name() << "\n";
                        break;
                    };
                };
            };
                
        } else {
                
            // Not binary
            
            // Go through all units
            for (auto unit=0; unit<_box_length; unit++) {
                
                // Write pos
                f << _look_up_pos(unit) << " ";

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

	void Lattice1DFullyVisible::read_from_file(MCType chain, int i_chain, std::string fname, bool binary)
	{
		// Clear the layer
        set_empty_all_units(chain, i_chain);

		// Open
		std::ifstream f;
		f.open(fname);
		if (!f.is_open()) { // make sure we found it
			std::cerr << ">>> Error: Lattice1DFullyVisible::read_layer_from_file <<< could not find file: " << fname << std::endl;
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

        if (binary) {
            while (getline(f,line)) {
                if (line == "") { continue; };
                iss = std::istringstream(line);
                iss >> x >> sp;
                s = _look_up_unit(atoi(x.c_str()));
                occs_to_write_binary[sp].push_back(s);
                // Reset
                sp=""; x="";
            };
        } else if (!binary) {
            while (getline(f,line)) {
                if (line == "") { continue; };
                iss = std::istringstream(line);
                iss >> x;
                s = _look_up_unit(atoi(x.c_str()));
                for (auto i=0; i<_species_possible_vec.size(); i++) {
                    iss >> sp >> prob;
                    prob_val = atof(prob.c_str());
                    occs_to_write_prob[sp].push_back(std::make_pair(s,prob_val));
                    // Reset
                    sp=""; prob="";
                };
                // Reset
                x="";
            };
        };

        if (binary) {
            // Binary mode
            for (auto const &pr: occs_to_write_binary) {
                // find the species
                auto it = _species_possible_map.find(pr.first);
                if (it == _species_possible_map.end()) {
                    std::cerr << ">>> Lattice1DFullyVisible::read_from_file <<< Error: could not find species (binary): " << pr.first << std::endl;
                    exit(EXIT_FAILURE);
                };
                for (auto const &site: pr.second) {
                    _mc_chains[chain][i_chain][it->second][site] = 1.0;
                };
            };
        } else {
            // Prob mode
            for (auto const &pr: occs_to_write_prob) {
                // find the species
                auto it = _species_possible_map.find(pr.first);
                if (it == _species_possible_map.end()) {
                    std::cerr << ">>> Lattice1DFullyVisible::read_from_file <<< Error: could not find species (prob): " << pr.first << std::endl;
                    exit(EXIT_FAILURE);
                };
                for (auto const &site_pr: pr.second) {
                    _mc_chains[chain][i_chain][it->second][site_pr.first] = site_pr.second;
                };
            };
        };
        
		f.close();
	};

    // ***************
    // MARK: - Activate layer steps
    // ***************
    
    // Calculate activation given layer above or below
    void Lattice1DFullyVisible::activate_calculate(MCType chain, int idx_start, int skip) {
        
        double bias;
        for (auto i_chain=0; i_chain<_no_markov_chains.at(chain); i_chain++) {
            for (auto sp: _species_possible_vec) {
                bias = get_bias(sp);
                for (auto x=idx_start; x<_box_length; x+=skip) {
                    // bias term
                    _mc_chains_act[chain][i_chain][sp][x] = bias;
                    
                    // ixns
                    if (_conn_2) {
                        // Backward
                        if (x != 0) {
                            for (auto &sp_other: _species_possible_vec) {
                                _mc_chains_act[chain][i_chain][sp][x] += get_ixn_2(sp, sp_other) * _mc_chains.at(chain).at(i_chain).at(sp_other).at(x-1);
                            };
                        };
                        // Forward
                        if (x != _box_length-1) {
                            for (auto &sp_other: _species_possible_vec) {
                                _mc_chains_act[chain][i_chain][sp][x] += get_ixn_2(sp, sp_other) * _mc_chains.at(chain).at(i_chain).at(sp_other).at(x+1);
                            };
                        };
                    };
                    
                    if (_conn_3) {
                        // I am third [O O X]
                        if (x != 0 && x != 1) {
                            for (auto &sp_other_1: _species_possible_vec) {
                                for (auto &sp_other_2: _species_possible_vec) {
                                    _mc_chains_act[chain][i_chain][sp][x] += get_ixn_3(sp_other_1,sp_other_2,sp) * _mc_chains.at(chain).at(i_chain).at(sp_other_1).at(x-2) * _mc_chains.at(chain).at(i_chain).at(sp_other_2).at(x-1);
                                };
                            };
                        };
                        // I am middle [O X O]
                        if (x != 0 && x != _box_length-1) {
                            for (auto &sp_other_1: _species_possible_vec) {
                                for (auto &sp_other_3: _species_possible_vec) {
                                    _mc_chains_act[chain][i_chain][sp][x] += get_ixn_3(sp_other_1,sp,sp_other_3) * _mc_chains.at(chain).at(i_chain).at(sp_other_1).at(x-1) * _mc_chains.at(chain).at(i_chain).at(sp_other_3).at(x+1);
                                };
                            };
                        };
                        // I am first [X O O]
                        if (x != _box_length-2 && x != _box_length-1) {
                            for (auto &sp_other_2: _species_possible_vec) {
                                for (auto &sp_other_3: _species_possible_vec) {
                                    _mc_chains_act[chain][i_chain][sp][x] += get_ixn_3(sp,sp_other_2,sp_other_3) * _mc_chains.at(chain).at(i_chain).at(sp_other_2).at(x+1) * _mc_chains.at(chain).at(i_chain).at(sp_other_3).at(x+2);
                                };
                            };
                        };
                    };

                };
            };
        };
    };

    // (2) Convert activations to probs
    void Lattice1DFullyVisible::activate_convert_to_probs(MCType chain, bool binary, int idx_start, int skip) {
        
        // Starts at 1.0 = exp(0) for empty
        _pst_prop.fill(arma::fill::ones);

        // Units
        for (auto x=idx_start; x<_box_length; x+=skip) {

            // All chains
            for (auto i_chain=0; i_chain<_no_markov_chains.at(chain); i_chain++) {
        
                // Convert activations to propensities via exp
                // Also calculate total propensity
                _pst_prop.fill(arma::fill::ones);
                for (auto sp: _species_possible_vec) {
                    _mc_chains_act[chain][i_chain][sp][x] = exp(_mc_chains_act[chain][i_chain][sp][x]);
                    _pst_prop[x] += _mc_chains_act[chain][i_chain][sp][x];
                };
                
                // Divide by total to normalize
                for (auto sp: _species_possible_vec) {
                    _mc_chains_act[chain][i_chain][sp][x] /= _pst_prop[x];
                };
            };
            
            // Sample if binary
            if (binary) {
                for (auto i_chain=0; i_chain<_no_markov_chains.at(chain); i_chain++) {
                    binarize_unit(chain, i_chain, x, true);
                };
            };
        };
    };
    
    // (3) Commit the new probabilities
    void Lattice1DFullyVisible::activate_committ(MCType chain, int idx_start, int skip) {
        // Units
        for (auto x=idx_start; x<_box_length; x+=skip) {
            // All chains
            for (auto i_chain=0; i_chain<_no_markov_chains.at(chain); i_chain++) {
                for (auto sp: _species_possible_vec) {
                    _mc_chains[chain][i_chain][sp][x] = _mc_chains_act[chain][i_chain][sp][x];
                };
            };
        };
    };
    
    // ***************
    // MARK: - Mean field / gibbs sampling
    // ***************
    
    // Sample
    void Lattice1DFullyVisible::gibbs_sampling_step() {
        
        int idx_start_min = 0;
        int idx_start_max = 0;
        int skip = 1;
        
        if (_conn_2) {
            idx_start_max = 1;
            skip = 2;
        };
        if (_conn_3) {
            idx_start_max = 2;
            skip = 3;
        };

        bool binary=true;
        for (int idx_start=idx_start_min; idx_start<=idx_start_max; idx_start++) {
            activate_calculate(MCType::ASLEEP, idx_start, skip);
            activate_convert_to_probs(MCType::ASLEEP, binary, idx_start, skip);
            activate_committ(MCType::ASLEEP, idx_start, skip);
        };
    };
    
    // ***************
    // MARK: - Reap moments
    // ***************
    
    void Lattice1DFullyVisible::reap_ixn_moment_diffs() {
        
        // O2 ixns
        Sptr sp1, sp2, sp3;
        std::shared_ptr<MomentDiff> moment;
        if (_conn_2) {
            for (auto &o2a: _o2_ixn_dict) {
                sp1 = o2a.first;
                for (auto &o2b: o2a.second) {
                    sp2 = o2b.first;
                    moment = o2b.second->get_moment_diff();
                    
                    // Reset moments
                    if (!moment->get_is_awake_moment_fixed()) {
                        moment->reset_moment(MCType::AWAKE);
                    };
                    moment->reset_moment(MCType::ASLEEP);

                    for (auto x=0; x<_box_length-1; x++) {
                        if (!moment->get_is_awake_moment_fixed()) {
                            for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::AWAKE); i_chain++) {
                                moment->increment_moment(MCType::AWAKE, _mc_chains.at(MCType::AWAKE).at(i_chain).at(sp1).at(x) * _mc_chains.at(MCType::AWAKE).at(i_chain).at(sp2).at(x+1) / _no_markov_chains.at(MCType::AWAKE));
                                // Note: the other order is automatically covered?
                            };
                        };
                        
                        for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::ASLEEP); i_chain++) {
                            moment->increment_moment(MCType::ASLEEP, _mc_chains.at(MCType::ASLEEP).at(i_chain).at(sp1).at(x) * _mc_chains.at(MCType::ASLEEP).at(i_chain).at(sp2).at(x+1) / _no_markov_chains.at(MCType::ASLEEP));
                            // Note: the other order is automatically covered?
                        };
                    };
                };
            };
        };
        
        // O3 ixns
        if (_conn_3) {
            for (auto &o3a: _o3_ixn_dict) {
                sp1 = o3a.first;
                for (auto &o3b: o3a.second) {
                    sp2 = o3b.first;
                    for (auto &o3c: o3b.second) {
                        sp3 = o3c.first;
                        moment = o3c.second->get_moment_diff();
                        
                        // Reset moments
                        if (!moment->get_is_awake_moment_fixed()) {
                            moment->reset_moment(MCType::AWAKE);
                        };
                        moment->reset_moment(MCType::ASLEEP);
                        
                        for (auto x=0; x<_box_length-2; x++) {
                            if (!moment->get_is_awake_moment_fixed()) {
                                for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::AWAKE); i_chain++) {
                                    moment->increment_moment(MCType::AWAKE, _mc_chains.at(MCType::AWAKE).at(i_chain).at(sp1).at(x) * _mc_chains.at(MCType::AWAKE).at(i_chain).at(sp2).at(x+1) * _mc_chains.at(MCType::AWAKE).at(i_chain).at(sp3).at(x+2) / _no_markov_chains.at(MCType::AWAKE));
                                    // Note: the other order is automatically covered?
                                };
                            };
                            
                            for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::ASLEEP); i_chain++) {
                                moment->increment_moment(MCType::ASLEEP, _mc_chains.at(MCType::ASLEEP).at(i_chain).at(sp1).at(x) * _mc_chains.at(MCType::ASLEEP).at(i_chain).at(sp2).at(x+1) * _mc_chains.at(MCType::ASLEEP).at(i_chain).at(sp3).at(x+2) / _no_markov_chains.at(MCType::ASLEEP));
                                // Note: the other order is automatically covered?
                            };
                        };
                    };
                };
            };
        };

            
        // Reap biases
        for (auto &bsp: _bias_dict) {
            sp1 = bsp.first;
            moment = bsp.second->get_moment_diff();
            
            if (!moment->get_is_awake_moment_fixed()) {
                moment->reset_moment(MCType::AWAKE);
            };
            moment->reset_moment(MCType::ASLEEP);

            if (!moment->get_is_awake_moment_fixed()) {
                for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::AWAKE); i_chain++) {
                    
                    moment->increment_moment(MCType::AWAKE, arma::accu(_mc_chains.at(MCType::AWAKE).at(i_chain).at(sp1)) / _no_markov_chains.at(MCType::AWAKE));
                };
            };
            
            for (auto i_chain=0; i_chain<_no_markov_chains.at(MCType::ASLEEP); i_chain++) {
                moment->increment_moment(MCType::ASLEEP, arma::accu(_mc_chains.at(MCType::ASLEEP).at(i_chain).at(sp1)) / _no_markov_chains.at(MCType::ASLEEP));
            };
        };
    };
    
    // ***************
    // MARK: - Add ixn to all ixns vec
    // ***************
    
    // Add ixn to all ixns vec
    void Lattice1DFullyVisible::_add_to_all_ixns_vec(Iptr ixn) {
    
        auto it = std::find(_all_ixns.begin(),_all_ixns.end(),ixn);
        if (it == _all_ixns.end()) {
            _all_ixns.push_back(ixn);
        } else {
            /*
            std::cerr << ">>> Lattice1DFullyVisible::_add_to_all_ixns_vec <<< reusing ixns currently not supported because it is not treated correctly in the reap function!" << std::endl;
            exit(EXIT_FAILURE);
             */
        };
    };
    
    // ***************
    // MARK: - Wake/sleep loop
    // ***************
    
    void Lattice1DFullyVisible::wake_sleep_loop_cd(int i_opt_step, int no_cd_steps, std::vector<FName> &fnames, OptionsWakeSleep_1DFV_CD options) {
        
        // AWAKE PHASE
        
        clock_t t0 = clock();
        
        // Read in the batch
        for (int i_chain=0; i_chain<_no_markov_chains[MCType::AWAKE]; i_chain++)
        {
            read_from_file(MCType::AWAKE, i_chain, fnames[i_chain].name, fnames[i_chain].binary);
        };
        _mc_chains[MCType::ASLEEP] = _mc_chains.at(MCType::AWAKE);
        
        clock_t t1 = clock();
        
        // Run CD sampling
        
        // Sample vis, hidden repeadedly
        for (int i_sampling_step=0; i_sampling_step<no_cd_steps; i_sampling_step++)
        {
            gibbs_sampling_step();
        };
        
        clock_t t2 = clock();

        if (options.verbose_timing) {
            double dt1 = (t1-t0)  / (double) CLOCKS_PER_SEC;
            double dt2 = (t2-t1)  / (double) CLOCKS_PER_SEC;
            double dt_tot = dt1 + dt2;
            std::cout << "[time " << dt_tot << "] [read " << dt1/dt_tot << "] [gibbs " << dt2/dt_tot << "]" << std::endl;
        };
    };
    
    // ***************
    // MARK: - Counts
    // ***************
    
    double Lattice1DFullyVisible::get_count(MCType chain, int i_chain, Sptr sp) const {
        return arma::accu(_mc_chains.at(chain).at(i_chain).at(sp));
    };
    double Lattice1DFullyVisible::get_count(MCType chain, int i_chain, Sptr sp1, Sptr sp2) const {
        double count=0.;
        for (auto x=0; x<_box_length-1; x++) {
            count += _mc_chains.at(chain).at(i_chain).at(sp1).at(x) * _mc_chains.at(chain).at(i_chain).at(sp2).at(x+1);
        };
        if (sp1 != sp2) {
            for (auto x=0; x<_box_length-1; x++) {
                count += _mc_chains.at(chain).at(i_chain).at(sp2).at(x) * _mc_chains.at(chain).at(i_chain).at(sp1).at(x+1);
            };
        };
        return count;
    };
    double Lattice1DFullyVisible::get_count(MCType chain, int i_chain, Sptr sp1, Sptr sp2, Sptr sp3) const {
        double count=0.;
        for (auto x=0; x<_box_length-2; x++) {
            count += _mc_chains.at(chain).at(i_chain).at(sp1).at(x) * _mc_chains.at(chain).at(i_chain).at(sp2).at(x+1) * _mc_chains.at(chain).at(i_chain).at(sp3).at(x+2);
        };
        if (!((sp1 == sp2) && (sp2 == sp3))) {
            for (auto x=0; x<_box_length-2; x++) {
                count += _mc_chains.at(chain).at(i_chain).at(sp3).at(x) * _mc_chains.at(chain).at(i_chain).at(sp2).at(x+1) * _mc_chains.at(chain).at(i_chain).at(sp1).at(x+2);
            };
        };
        return count;
    };
    double Lattice1DFullyVisible::get_count(MCType chain, int i_chain, Sptr sp1, Sptr sp2, Sptr sp3, Sptr sp4) const {
        double count=0.;
        for (auto x=0; x<_box_length-3; x++) {
            count += _mc_chains.at(chain).at(i_chain).at(sp1).at(x) * _mc_chains.at(chain).at(i_chain).at(sp2).at(x+1) * _mc_chains.at(chain).at(i_chain).at(sp3).at(x+2) * _mc_chains.at(chain).at(i_chain).at(sp4).at(x+3);
        };
        if (!((sp1 == sp2) && (sp2 == sp3) && (sp3 == sp4))) {
            for (auto x=0; x<_box_length-3; x++) {
                count += _mc_chains.at(chain).at(i_chain).at(sp4).at(x) * _mc_chains.at(chain).at(i_chain).at(sp3).at(x+1) * _mc_chains.at(chain).at(i_chain).at(sp2).at(x+2) * _mc_chains.at(chain).at(i_chain).at(sp1).at(x+3);
            };
        };
        return count;
    };
};
