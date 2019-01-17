#include "../include/bmla_bits/opt_problem.hpp"

// Other headers
#include "../include/bmla_bits/ixn_param.hpp"
#include "../include/bmla_bits/lattice.hpp"
#include "../include/bmla_bits/moment.hpp"

#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>

/************************************
 * Namespace for bmla
 ************************************/

namespace bmla {
    
    /****************************************
     OptProblem
     ****************************************/
    
    /********************
     Constructor
     ********************/
    
    OptProblem::OptProblem(std::shared_ptr<Lattice> latt, std::vector<std::shared_ptr<IxnParam>> ixn_params) {
        _latt = latt;
        _ixn_params = ixn_params;
    };
    OptProblem::OptProblem(const OptProblem& other) {
        _copy(other);
    };
    OptProblem::OptProblem(OptProblem&& other) {
        _move(other);
    };
    OptProblem& OptProblem::operator=(const OptProblem& other) {
        if (this != &other) {
            _clean_up();
            _copy(other);
        };
        return *this;
    };
    OptProblem& OptProblem::operator=(OptProblem&& other) {
        if (this != &other) {
            _clean_up();
            _move(other);
        };
        return *this;
    };
    OptProblem::~OptProblem() {
        _clean_up();
    };
    
    void OptProblem::_clean_up() {};
    void OptProblem::_move(OptProblem &other) {
        _ixn_params = std::move(other._ixn_params);
        _latt = std::move(other._latt);
    };
    void OptProblem::_copy(const OptProblem& other) {
        _ixn_params = other._ixn_params;
        if (other._latt) {
            _latt = std::make_shared<Lattice>(*other._latt);
        };
    };
    
    /********************
     Init structures
     ********************/
    
    void OptProblem::init_structures(int batch_size, int no_markov_chains) {
        
        // Ixn params
        for (auto &ixn_param: _ixn_params) {
            
            // Moment
            ixn_param->get_moment()->set_batch_size(batch_size);
            ixn_param->get_moment()->set_no_markov_chains(no_markov_chains);
        };
        
        // Lattice
        _latt->set_no_markov_chains(no_markov_chains);
    };
};
