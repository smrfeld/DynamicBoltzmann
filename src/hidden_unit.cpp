#include "math.h"
#include <iostream>
#include "../include/general.hpp"
#include "ixn_param_traj.hpp" // includes hidden unit hpp

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
		_bias = nullptr;
		_t_opt_ptr = nullptr;
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
		_bias = nullptr;
		_t_opt_ptr = nullptr;
	};
	void HiddenUnit::_copy(const HiddenUnit& other) {
		_n_conn = other._n_conn;
		_conn = other._conn;
		_sp = other._sp;
		_val = other._val;
		_bias = other._bias;
		_t_opt_ptr = other._t_opt_ptr;
	};

	/********************
	Getters
	********************/

	double HiddenUnit::get() const {
		return _val;
	};

	/********************
	Set the time ptr
	********************/

	void HiddenUnit::set_t_opt_ptr(int *t_opt_ptr) {
		_t_opt_ptr = t_opt_ptr;
	};

	/********************
	Set the bias
	********************/

	void HiddenUnit::set_bias(IxnParamTraj *ip) {
		_bias = ip;
	};

	/********************
	Activate
	********************/

	void HiddenUnit::activate(bool binary) {
		double act = 0.0;

		// Bias
		if (_bias) {
			act += _bias->get_at_time(*_t_opt_ptr);
		};

		// Go through all connected neurons
		for (auto c: _conn) {
			if (c->sp == _sp) { // Check that this is the species that I love
				// act += c->sp->w(); // Gets the weight of this connection
				// std::cout << "Activating hidden:   " << c->sp->w() << std::endl;
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
		return 1.0 / (1.0 + exp(-x));
	};	
};