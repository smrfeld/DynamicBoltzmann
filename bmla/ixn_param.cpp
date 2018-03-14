#include "ixn_param.hpp"
#include "../general.hpp"
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
		_val = other._val;
		_val_guess = other._val_guess;
		_asleep = other._asleep;
		_awake = other._awake;
	};
	void IxnParam::_copy(IxnParam&& other)
	{
		_name = other._name;
		_type = other._type;
		_sp1 = other._sp1;
		_sp2 = other._sp2;
		_val = other._val;
		_val_guess = other._val_guess;
		_asleep = other._asleep;
		_awake = other._awake;
		// Clear other
		other._name = "";
		other._sp1 = nullptr;
		other._sp2 = nullptr;
		other._val = 0.;
		other._val_guess = 0.;
		other._asleep = 0.;
		other._awake = 0.;
	};
	void IxnParam::_clean_up() {};

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
			if (_sp1) {
				if (moment_type==AWAKE) {
					_awake += 1. * _sp1->count() / batch_size;
				} else if (moment_type==ASLEEP) {
					_asleep += 1. * _sp1->count() / batch_size;
				};
			};
		} else if (_type == Jp) {
			if (_sp1 && _sp2) {
				if (moment_type==AWAKE) {
					_awake += 1. * _sp1->nn_count(_sp2) / batch_size;
				} else if (moment_type==ASLEEP) {
					_asleep += 1. * _sp1->nn_count(_sp2) / batch_size;
				};
			};
		};	
	};

	double IxnParam::moments_diff() const {
		return (_awake - _asleep);
	};

};





