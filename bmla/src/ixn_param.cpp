#include "ixn_param.hpp"
#include "../include/general.hpp"
#include <iostream>
#include <fstream>
#include "math.h"
#include <algorithm>

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

	IxnParam::IxnParam(std::string name, IxnParamType type, Species *sp, double val_guess) : IxnParam(name,type,sp,nullptr,nullptr,val_guess) {};
	IxnParam::IxnParam(std::string name, IxnParamType type, Species *sp1, Species *sp2, double val_guess) : IxnParam(name,type,sp1,sp2,nullptr,val_guess) {};
	IxnParam::IxnParam(std::string name, IxnParamType type, Species *sp1, Species *sp2, Species *sp3, double val_guess)	
	{
		_name = name;
		_type = type;
		_sp1 = sp1;
		_sp2 = sp2;
		_sp3 = sp3;
		_val_guess = val_guess;

		_val = val_guess;

		_track_soln_traj = false;

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
		_sp3 = other._sp3;
		_conns = other._conns;
		_hidden_units = other._hidden_units;
		_val = other._val;
		_val_guess = other._val_guess;
		_track_soln_traj = other._track_soln_traj;
		_soln_traj = other._soln_traj;
		_asleep = other._asleep;
		_awake = other._awake;
	};
	void IxnParam::_reset()
	{
		_name = "";
		_sp1 = nullptr;
		_sp2 = nullptr;
		_sp3 = nullptr;
		_conns.clear();
		_hidden_units.clear();
		_val = 0.;
		_val_guess = 0.;
		_track_soln_traj = false;
		_soln_traj.clear();
		_asleep = 0.;
		_awake = 0.;
	};
	void IxnParam::_clean_up() {};

	/********************
	Check if this ixn param is...
	********************/

	bool IxnParam::is_type_with_species(IxnParamType type, std::string s) const {
		if (_type == type) {
			if (_sp1->name() == s) {
				return true;
			};
		};
		return false;
	};
	bool IxnParam::is_type_with_species(IxnParamType type, std::string s1, std::string s2) const {
		if (_type == type) {
			if ((_sp1->name() == s1 && _sp2->name() == s2) || (_sp1->name() == s2 && _sp2->name() == s1)) {
				return true;
			};
		};
		return false;
	};
	bool IxnParam::is_type_with_species(IxnParamType type, std::string s1, std::string s2, std::string s3) const{
		if (_type == type) {
			if ((_sp1->name() == s1 && _sp2->name() == s2 && _sp3->name() == s3) 
				|| (_sp1->name() == s3 && _sp2->name() == s1 && _sp3->name() == s2)
				|| (_sp1->name() == s2 && _sp2->name() == s3 && _sp3->name() == s1)
				) {
				return true;
			};
		};
		return false;
	};

	/********************
	Add a visible->hidden unit connection
	If we are Wp
	********************/

	void IxnParam::add_visible_hidden_connection(Site *sptr, HiddenUnit *hup) {
		if (_type != Wp) {
			std::cerr << "ERROR: adding visible to hidden connection to wrong type of IxnParam" << std::endl;
			exit(EXIT_FAILURE);
		};
		_conns.push_back(std::make_pair(sptr,hup));
	};

	/********************
	Add a hidden unit to monitor
	If we are Bp
	********************/

	void IxnParam::add_hidden_unit(HiddenUnit *hup) {
		if (_type != Bp) {
			std::cerr << "ERROR: adding hidden unit to wrong type of IxnParam" << std::endl;
			exit(EXIT_FAILURE);
		};
		_hidden_units.push_back(hup);
	};

	/********************
	Update
	********************/

	void IxnParam::update(double dopt, bool l2_reg, double lambda) {
		// Store the old if needed
		if (_track_soln_traj) {
			_soln_traj.push_back(_val);
		};

		// Get new
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
	void IxnParam::set_val(double val) {
		_val = val;
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
		} else if (_type == Kp) {
			if (moment_type==AWAKE) {
				_awake += 1. * _sp1->triplet_count(_sp2,_sp3) / batch_size;
			} else if (moment_type==ASLEEP) {
				_asleep += 1. * _sp1->triplet_count(_sp2,_sp3) / batch_size;
			};
		} else if (_type == Wp) {
			// Run through all connections
			double inc = 0.0;
			for (auto c: _conns) {
				/*
				std::cout << "visible: " << c.first->x << " " << c.first->y << " " << c.first->z << " hidden: ";
				c.second->print_conns(false);
				std::cout << " vis val: " << c.first->get_prob(_sp1) << " hidden val: " << c.second->get() << std::endl;
				*/
				
				inc += c.first->get_prob(_sp1) * c.second->get(); // v * h
			};
			if (moment_type==AWAKE) {
				_awake += 1. * inc / batch_size;
			} else if (moment_type==ASLEEP) {
				_asleep += 1. * inc / batch_size;
			};
		} else if (_type == Bp) {
			// Run through all hidden units that are monitored
			double inc = 0.0;
			for (auto hup: _hidden_units) {
				inc += hup->get();
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

	/********************
	Store soln traj
	********************/

	void IxnParam::set_track_soln_traj(bool flag) {
		_track_soln_traj = flag;
	};
	void IxnParam::reset_soln_traj() {
		_soln_traj.clear();
	};
	double IxnParam::get_ave() const {
		return get_ave(_soln_traj.size());
	};
	double IxnParam::get_ave(int last_n_steps) const {
		if (!_track_soln_traj) {
			std::cerr << "Error! Requesting average but solution was not tracked. Turn on track_soln_traj option!" << std::endl;
			exit(EXIT_FAILURE);
		};

		double ave=0.0;
		int i_start = std::max(int(_soln_traj.size()-last_n_steps),0);
		int n = _soln_traj.size() - i_start;
		for (auto i=i_start; i<_soln_traj.size(); i++) {
			ave += _soln_traj[i];
		};
		return ave / n;
	};
};





