#include "hidden_unit.hpp"
#include "math.h"
#include <iostream>
#include "../include/general.hpp"

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

	HiddenUnit::HiddenUnit(std::vector<Site*> conn, Species *sp)
	{
		_n_conn = conn.size();
		_conn = conn;
		_sp = sp;
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
		_n_conn = 0;
		_conn.clear();
		_sp = nullptr;
		_val = 0.;
	};
	void HiddenUnit::_copy(const HiddenUnit& other) {
		_n_conn = other._n_conn;
		_conn = other._conn;
		_sp = other._sp;
		_val = other._val;
	};

	/********************
	Print connections
	********************/

	void HiddenUnit::print_conns(bool newline) const {
		std::cout << "Connected to: ";
		for (auto sp: _conn) {
			std::cout << "(" << sp->x << "," << sp->y << "," << sp->z << ") ";
		};
		if (newline) {
			std::cout << std::endl;
		};
	};

	/********************
	Getters
	********************/

	double HiddenUnit::get() const {
		return _val;
	};

	/********************
	Activate
	********************/

	void HiddenUnit::activate(bool binary) {
		// Go through all connected neurons
		//std::cout << "Activating..." << std::endl;
		double act = 0.0;
		for (auto c: _conn) {
			if (c->binary) {
				// Binary
				if (c->sp == _sp) { // Check that this is the species that I love
					act += _sp->w(); // Gets the weight of this connection
				};
			} else {
				// Probabilistic
				act += _sp->w() * c->get_prob(_sp); // weight * visible value for this species
			};
		};

		// Pass through sigmoid
		_val = _sigma(act);
		if (binary) {
			// 1 if probability = _val
			if (randD(0.0,1.0) <= _val) {
				_val = 1.;
			} else {
				_val = 0.;
			};
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