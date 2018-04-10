#include "ixn_param.hpp"
#include "../include/general.hpp"
#include <iostream>
#include <fstream>
#include "math.h"

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Ixn Param
	****************************************/

	/********************
	Constructor
	********************/

	IxnParam::IxnParam(std::string name, IxnParamType type, Species *sp, double val_guess) : IxnParam(name,type,sp,nullptr,val_guess) {};
	IxnParam::IxnParam(std::string name, IxnParamType type, Species *sp1, Species *sp2, double val_guess)
	{
		_name = name;
		_type = type;
		_sp1 = sp1;
		_sp2 = sp2;
		_val_guess = val_guess;

		_val = val_guess;

		_asleep = 0.0;
		_awake = 0.0;
	};
	IxnParam::IxnParam(const IxnParam& other) {
		_copy(other);
	};
	IxnParam::IxnParam(IxnParam&& other) {
		_copy(other);
		other._reset();
	};
	IxnParam& IxnParam::operator=(const IxnParam& other) {
		if (this != &other)
		{
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	IxnParam& IxnParam::operator=(IxnParam&& other) {
		if (this != &other)
		{
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
	};
	IxnParam::~IxnParam() {
		_clean_up();
	};
	void IxnParam::_copy(const IxnParam& other)
	{
		_name = other._name;
		_type = other._type;
		_sp1 = other._sp1;
		_sp2 = other._sp2;
		_conns = other._conns;
		_val = other._val;
		_val_guess = other._val_guess;
		_asleep = other._asleep;
		_awake = other._awake;
	};
	void IxnParam::_reset()
	{
		_name = "";
		_sp1 = nullptr;
		_sp2 = nullptr;
		_conns.clear();
		_val = 0.;
		_val_guess = 0.;
		_asleep = 0.;
		_awake = 0.;
	};
	void IxnParam::_clean_up() {};

	/********************
	Check if this ixn param is...
	********************/

	bool IxnParam::is_w_with_species(std::string species_name) const {
		if (_type == Wp && _sp1->name() == species_name) {
			return true;
		} else {
			return false;
		};
	};

	bool IxnParam::is_j_with_species(std::string species_name_1, std::string species_name_2) const {
		if (_type == Jp) {
			if ((_sp1->name() == species_name_1 && _sp2->name() == species_name_2) || (_sp1->name() == species_name_2 && _sp2->name() == species_name_1)) {
				return true;
			};
		};
		return false;
	};

	/********************
	Add a visible->hidden unit connection
	********************/

	void IxnParam::add_visible_hidden_connection(Site *sptr, HiddenUnit *hup) {
		if (_type != Wp) {
			std::cerr << "ERROR: adding visible to hidden connection to wrong type of IxnParam" << std::endl;
			exit(EXIT_FAILURE);
		};
		_conns.push_back(std::make_pair(sptr,hup));
	};

	/********************
	Update
	********************/

	void IxnParam::update(double dopt, bool l2_reg, double lambda) {
		_val += dopt * moments_diff();
		if (l2_reg) {
			if (_val > 0) {
				_val -= lambda * 2. * abs(_val);
			} else {
				_val += lambda * 2. * abs(_val);
			};
		};
	};

	/********************
	Getters/setters
	********************/

	std::string IxnParam::name() const {
		return _name;
	};

	double IxnParam::get() const
	{
		return _val;
	};
	void IxnParam::set_guess(double guess) {
		_val_guess = guess;
	};

	void IxnParam::reset() {
		_val = _val_guess;
	};

	/********************
	Moments from lattice
	********************/

	double IxnParam::get_moment(MomentType moment_type) const {
		if (moment_type == AWAKE) {
			return _awake;
		} else {
			return _asleep;
		};
	};

	void IxnParam::moments_reset(MomentType moment_type) 
	{
		if (moment_type==AWAKE) {
			_awake = 0.;
		} else {
			_asleep = 0.;
		};
	};
	void IxnParam::moments_retrieve(MomentType moment_type) {
		moments_retrieve(moment_type, 1);
	};
	void IxnParam::moments_retrieve(MomentType moment_type, int batch_size)
	{
		if (_type == Hp) {
			if (moment_type==AWAKE) {
				_awake += 1. * _sp1->count() / batch_size;
			} else if (moment_type==ASLEEP) {
				_asleep += 1. * _sp1->count() / batch_size;
			};
		} else if (_type == Jp) {
			if (moment_type==AWAKE) {
				_awake += 1. * _sp1->nn_count(_sp2) / batch_size;
			} else if (moment_type==ASLEEP) {
				_asleep += 1. * _sp1->nn_count(_sp2) / batch_size;
			};
		} else if (_type == Wp) {
			double inc = 0.0;
			// Run through all connections
			for (auto c: _conns) {
				// Is it binary or not?
				if (c.first->binary) {
					// Binary
					/*
					std::cout << "visible: " << c.first->x << " " << c.first->y << " " << c.first->z << " hidden: ";
					c.second->print_conns(false);
					std::cout << " vis val: ";
					if (c.first->sp == _sp1) { // site occupied with the correct species
						std::cout << 1;
					} else {
						std::cout << 0;
					};
					std::cout << " hidden val: " << c.second->get() << std::endl;
					*/

					if (c.first->sp == _sp1) { // site occupied with the correct species
						inc += c.second->get(); // v * h
					};
				} else {
					// Probabilistic
					// Use the prob for the species we love
					/*
					std::cout << "visible: " << c.first->x << " " << c.first->y << " " << c.first->z << " hidden: ";
					c.second->print_conns(false);
					std::cout << " vis val: " << c.first->get_prob(_sp1) << " hidden val: " << c.second->get() << std::endl;
					*/
					
					inc += c.first->get_prob(_sp1) * c.second->get(); // v * h
				};
			};
			if (moment_type==AWAKE) {
				_awake += 1. * inc / batch_size;
			} else if (moment_type==ASLEEP) {
				_asleep += 1. * inc / batch_size;
			};
		};
	};

	double IxnParam::moments_diff() const {
		return (_awake - _asleep);
	};

};





