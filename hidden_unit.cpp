#include "hidden_unit.hpp"
#include "math.h"

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

	HiddenUnit::HiddenUnit() : HiddenUnit(std::vector<Site*>(), nullptr) {};
	HiddenUnit::HiddenUnit(std::vector<Site*> conn) : HiddenUnit(conn, nullptr) {};
	HiddenUnit::HiddenUnit(std::vector<Site*> conn, IxnParamTraj* weight)
	{
		_n_conn = conn.size();
		_conn = conn;
		_weight = weight;
		_t_opt_ptr = nullptr;
		_val = 0.0;
	};
	HiddenUnit::HiddenUnit(const HiddenUnit& other)
	{
		_copy(other);
	};
	HiddenUnit::HiddenUnit(HiddenUnit&& other)
	{
		_copy(other);
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
	void HiddenUnit::_copy(const HiddenUnit& other) {
		_n_conn = other._n_conn;
		_conn = other._conn;
		_weight = other._weight;
		_t_opt_ptr = other._t_opt_ptr;
		_val = other._val;
	};
	void HiddenUnit::_copy(HiddenUnit&& other) {
		_n_conn = other._n_conn;
		_conn = other._conn;
		_weight = other._weight;
		_t_opt_ptr = other._t_opt_ptr;
		_val = other._val;
		// Clear other
		other._n_conn = 0;
		other._conn.clear();
		other._weight = nullptr;
		other._t_opt_ptr = nullptr;
		other._val = 0.;
	};

	/********************
	Add connections
	********************/

	void HiddenUnit::add_conn(Site* site_ptr) {
		_conn.push_back(site_ptr);
	};

	/********************
	Set pointer to the opt time variable
	********************/

	void HiddenUnit::set_opt_time_ptr(int *t_opt_ptr) {
		_t_opt_ptr = t_opt_ptr;
	};

	/********************
	Activate
	********************/

	void HiddenUnit::activate(bool binary) {
		// Go through all connected neurons
		double act = 0.0;
		for (auto c: _conn) {
			if (c->sp != nullptr) {
				act += 1.0;
			};
		};
		// Multiply by weight
		act *= _weight->get_at_time(*_t_opt_ptr);
		// Pass through sigmoid
		_val = _sigma(act);
		if (binary) {
			if (_val >= 0.5) {
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