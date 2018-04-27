#include "math.h"
#include <iostream>
#include "../include/general.hpp"
#include "ixn_param.hpp" // includes hidden units

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

	HiddenUnit::HiddenUnit(std::vector<IxnParam*> bias)
	{
		_bias = bias;
		_val = 0.0;
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
		_conn.clear();
		_bias.clear();
		_val = 0.;
	};
	void HiddenUnit::_copy(const HiddenUnit& other) {
		_conn = other._conn;
		_bias = other._bias;
		_val = other._val;
	};

	/********************
	Add a connection
	********************/

	void HiddenUnit::add_connection(ConnectionVH* conn) {
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

	double HiddenUnit::get() const {
		return _val;
	};

	/********************
	Set the bias
	********************/

	void HiddenUnit::add_bias(IxnParam *ip) {
		_bias.push_back(ip);
	};

	/********************
	Activate
	********************/

	void HiddenUnit::activate(bool binary) {
		// Go through all connected neurons
		double act = 0.0;

		// Bias
		for (auto b: _bias) {
			act += b->get();
		};

		// Get act from conns
		for (auto c: _conn) {
			act += c->get_act_hidden();
		};
		/*
		// Weights - go through all connected visible units
		std::vector<Species*> sp_vec;
		for (auto c: _conn) {
			// Go through all ixn params associated with this connection
			for (auto ip: c.second) {
				// What species does this ixn param care about?
				sp_vec = ip->get_species(); // Most Wp only 1, or all species if all species apply
				for (auto sp: sp_vec) {
					// Weight * prob of this species at this site
					act += ip->get() * c.first->get_prob(sp);
				};
			};
		};
		*/
		
		// Pass through sigmoid -> probability
		_val = _sigma(act);

		// Binarize
		if (binary) {
			binarize();
		};
	};

	/********************
	Binarize
	********************/

	void HiddenUnit::binarize() {
		// 1 if probability = _val
		if (randD(0.0,1.0) <= _val) {
			_val = 1.;
		} else {
			_val = 0.;
		};
	};

	/********************
	PRIVATE - Activation function
	********************/

	double HiddenUnit::_sigma(double x) const {
		// Sigmoid
		return 1.0 / (1.0 + exp(-x));
	};	
};