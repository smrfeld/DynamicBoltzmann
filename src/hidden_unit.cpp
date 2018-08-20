#include "hidden_unit.hpp"

// Other headers
#include "../include/dynamicboltz_bits/general.hpp"
#include "lattice.hpp"
#include "species.hpp"
#include "ixn_param_traj.hpp"

#include "math.h"
#include <iostream>

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/************************************
	Hidden unit
	************************************/

	/********************
	Constructor
	********************/

	HiddenUnit::HiddenUnit()
	{
		_prob_empty = 1.0;
	};
	HiddenUnit::HiddenUnit(const HiddenUnit& other)
	{
		_copy(other);
	};
	HiddenUnit::HiddenUnit(HiddenUnit&& other)
	{
		_copy(other);
		other._reset();
	};
	HiddenUnit& HiddenUnit::operator=(const HiddenUnit& other)
	{
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	HiddenUnit& HiddenUnit::operator=(HiddenUnit&& other)
	{
		if (this != &other) {
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
	};
	HiddenUnit::~HiddenUnit()
	{
		_clean_up();
	};
	void HiddenUnit::_clean_up() {
		// Nothing....
	};
	void HiddenUnit::_reset() {
		_sp_possible.clear();
		_conn.clear();
		_bias.clear();
		_prob_empty = 1.0;
		_probs.clear();
	};
	void HiddenUnit::_copy(const HiddenUnit& other) {
		_sp_possible = other._sp_possible;
		_conn = other._conn;
		_bias = other._bias;
		_prob_empty = other._prob_empty;
		_probs = other._probs;
	};

	/********************
	Add possible species
	********************/

	void HiddenUnit::add_hidden_species_possibility(HiddenSpecies* sp) {
		_sp_possible.push_back(sp);
		_probs[sp] = 0.;
		_bias[sp] = std::vector<IxnParamTraj*>();
	};

	/********************
	Add a connection
	********************/

	void HiddenUnit::add_visible_hidden_conn(ConnectionVH* conn) {
		_conn.push_back(conn);
	};

	/********************
	Print connections
	********************/

	void HiddenUnit::print_conns(bool newline) const {
		/*
		std::cout << "Hidden unit connected to: ";
		for (auto c: _conn) {
			c->print();
			std::cout << "(" << c.first->x << "," << c.first->y << "," << c.first->z << ") with params: ";
			for (auto ip: c.second) {
				std::cout << ip->name() << " ";
			};
		};
		std::cout << " and biases: ";
		for (auto ip: _bias) {
			std::cout << ip->name() << " ";
		};
		if (newline) {
			std::cout << std::endl;
		};
		*/
	};

	/********************
	Getters
	********************/

	double HiddenUnit::get_prob(HiddenSpecies *hsp) const {
		if (!hsp) {
			// Prob empty
			return _prob_empty;
		};

		auto it = _probs.find(hsp);
		if (it != _probs.end()) {
			return it->second;
		} else {
			return 0.;
		};
	};

	/********************
	Set the bias
	********************/

	void HiddenUnit::add_bias(IxnParamTraj *ip) {
		// Get species involved
		std::vector<HiddenSpecies*> hsp_vec = ip->get_species_hidden_bias();

		// Go through species
		for (auto hsp: hsp_vec) {
			// Check this is a possible species for this site
			auto it = std::find(_sp_possible.begin(),_sp_possible.end(),hsp);
			if (it != _sp_possible.end()) {
				// Add
				_bias[hsp].push_back(ip);
			};
		};
	};

	/********************
	Activate
	********************/

	void HiddenUnit::activate(bool binary) {

		// Propensity vector
		std::vector<double> probs;
		std::vector<double> props;
		props.push_back(0.0);

		// Prop of empty
		probs.push_back(1.0);
		props.push_back(1.0);

		// Go through all possible species
		double energy;
		HiddenSpecies *hsp;
		for (auto hsp: _sp_possible) {

			// Reset energy
			energy = 0.0;

			// Get bias
			for (auto b: _bias[hsp]) {
				energy += b->get();
			};

			// Connections
			for (auto c: _conn) {
				energy += c->get_act_hidden(hsp);
			};

			// Store prop
			energy = exp(energy);
			props.push_back(props.back()+energy);
			probs.push_back(energy);
		};

		// Clear current
		for (auto hsp: _sp_possible) {
			_probs[hsp] = 0.0;
		};
		_prob_empty = 0.0;

		/*
		for (auto x: props) {
			std::cout << x << " ";
		};
		std::cout << "" << std::endl;
		*/
		
		// Store/sample
		int i_chosen;
		if (binary) {

			// Sample RV
			i_chosen = sample_prop_vec(props);

			if (i_chosen==0) {

				// Empty!
				_prob_empty = 1.0;

			} else {

				// Make species
				_probs[_sp_possible[i_chosen-1]] = 1.0;

			};

		} else {

			// Normalize probs
			double tot=0.0;
			for (auto pr: probs) {
				tot += pr;
			};

			// Write into species
			_prob_empty = probs[0]/tot;
			for (int i=0; i<_sp_possible.size(); i++) {
				_probs[_sp_possible[i]] = probs[i+1]/tot;
			};

		};
	};

	/********************
	Binarize
	********************/

	void HiddenUnit::binarize() {
		// 1 if probability = _val
		if (randD(0.0,1.0) <= _probs.begin()->second) {
			_probs[_probs.begin()->first] = 1.;
		} else {
			_probs[_probs.begin()->first] = 0.;
		};
	};
};