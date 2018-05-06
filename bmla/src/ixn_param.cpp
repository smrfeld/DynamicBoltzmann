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

	IxnParam::IxnParam(std::string name, IxnParamType type, double val_guess)	
	{
		_name = name;
		_type = type;

		_val_guess = val_guess;
		_val = val_guess;
		_val_nesterov = val_guess;

		_track_soln_traj = false;

		_asleep = 0.0;
		_awake = 0.0;

		_is_awake_fixed = false;
		_awake_fixed = 0.0;
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

		_sp_bias_visible = other._sp_bias_visible;
		_sp_bias_hidden = other._sp_bias_hidden;
		_sp_doublet = other._sp_doublet;
		_sp_triplet = other._sp_triplet;
		_sp_quartic = other._sp_quartic;
		_sp_conn = other._sp_conn;

		_conns = other._conns;
		_hidden_units = other._hidden_units;

		_val = other._val;
		_val_guess = other._val_guess;
		_val_nesterov = other._val_nesterov;
		_track_soln_traj = other._track_soln_traj;
		_soln_traj = other._soln_traj;
		_asleep = other._asleep;
		_awake = other._awake;
		_is_awake_fixed = other._is_awake_fixed;
		_awake_fixed = other._awake_fixed;
	};
	void IxnParam::_reset()
	{
		_name = "";
		
		_sp_bias_visible.clear();
		_sp_bias_hidden.clear();
		_sp_doublet.clear();
		_sp_triplet.clear();
		_sp_quartic.clear();
		_sp_conn.clear();

		_conns.clear();
		_hidden_units.clear();

		_val = 0.;
		_val_guess = 0.;
		_val_nesterov = 0.;
		_track_soln_traj = false;
		_soln_traj.clear();
		_asleep = 0.;
		_awake = 0.;

		_is_awake_fixed = false;
		_awake_fixed = 0.0;
	};
	void IxnParam::_clean_up() {};

	/********************
	Nesterov
	********************/

	// Move to the nesterov intermediate point
	void IxnParam::nesterov_move_to_intermediate_pt(int opt_step) {
		// Store the curren
		double curr = _val;
		// Move to the intermediate
		_val += (opt_step - 1.0) / (opt_step + 2.0) * (curr - _val_nesterov);
		// The current point will be the next old
		_val_nesterov = curr;
	};

	// Set prev nesterov
	void IxnParam::nesterov_set_prev_equal_curr() {
		_val_nesterov = _val;
	};


	/********************
	Add species
	********************/

	void IxnParam::add_species(Species *sp) {
		_sp_bias_visible.push_back(sp);
	};
	void IxnParam::add_species(HiddenSpecies *hsp) {
		_sp_bias_hidden.push_back(hsp);
	};
	void IxnParam::add_species(Species *sp1, Species *sp2) {
		_sp_doublet.push_back(Species2(sp1,sp2));
	};
	void IxnParam::add_species(Species *sp1, Species *sp2, Species *sp3) {
		_sp_triplet.push_back(Species3(sp1,sp2,sp3));
	};
	void IxnParam::add_species(Species *sp1, Species *sp2, Species *sp3, Species *sp4) {
		_sp_quartic.push_back(Species4(sp1,sp2,sp3,sp4));
	};
	void IxnParam::add_species(Species *sp, HiddenSpecies *hsp) {
		_sp_conn.push_back(SpeciesVH(sp,hsp));
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
	void IxnParam::set_fixed_awake_moment(double val) {
		_is_awake_fixed = true;
		_awake_fixed = val;
	};

	void IxnParam::reset_to_guess() {
		_val = _val_guess;
	};

	/********************
	Get species involved
	********************/

	std::vector<Species*> IxnParam::get_species_bias() const {
		if (_type != Hp) {
			std::cerr << "Error! Not a visible bias." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _sp_bias_visible;
	};
	std::vector<HiddenSpecies*> IxnParam::get_species_hidden_bias() const {
		if (_type != Bp) {
			std::cerr << "Error! Not a hidden bias." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _sp_bias_hidden;
	};
	std::vector<Species2> IxnParam::get_species_doublet() const {
		if (_type != Jp) {
			std::cerr << "Error! Not a J." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _sp_doublet;
	};
	std::vector<Species3> IxnParam::get_species_triplet() const {
		if (_type != Kp) {
			std::cerr << "Error! Not a K." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _sp_triplet;
	};
	std::vector<Species4> IxnParam::get_species_quartic() const {
		if (_type != Qp) {
			std::cerr << "Error! Not a Q." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _sp_quartic;
	};
	std::vector<SpeciesVH> IxnParam::get_species_conn() const {
		if (_type != Wp) {
			std::cerr << "Error! Not a conn W." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _sp_conn;
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
		if (_is_awake_fixed) {
			_awake = _awake_fixed;
			return;
		};
		if (_type == Hp) {
			for (auto sp: _sp_bias_visible) {
				if (moment_type==AWAKE) {
					_awake += 1. * sp->count() / batch_size;
				} else if (moment_type==ASLEEP) {
					_asleep += 1. * sp->count() / batch_size;
				};
			};
		} else if (_type == Jp) {
			for (auto spd: _sp_doublet) {
				if (moment_type==AWAKE) {
					_awake += 1. * spd.sp1->nn_count(spd.sp2) / batch_size;
				} else if (moment_type==ASLEEP) {
					_asleep += 1. * spd.sp1->nn_count(spd.sp2) / batch_size;
				};
			};
		} else if (_type == Kp) {
			for (auto spt: _sp_triplet) {
				if (moment_type==AWAKE) {
					_awake += 1. * spt.sp1->triplet_count(spt.sp2,spt.sp3) / batch_size;
				} else if (moment_type==ASLEEP) {
					_asleep += 1. * spt.sp1->triplet_count(spt.sp2,spt.sp3) / batch_size;
				};
			};
		} else if (_type == Qp) {
			for (auto spt: _sp_quartic) {
				if (moment_type==AWAKE) {
					_awake += 1. * spt.sp1->quartic_count(spt.sp2,spt.sp3,spt.sp4) / batch_size;
				} else if (moment_type==ASLEEP) {
					_asleep += 1. * spt.sp1->quartic_count(spt.sp2,spt.sp3,spt.sp4) / batch_size;
				};
			};
		} else if (_type == Wp) {
			// Run through all connections
			double inc = 0.0;
			for (auto c: _conns) {
				for (auto spvh: _sp_conn) {
					// v * h
					inc += c.first->get_prob(spvh.sp_visible) * c.second->get_prob(spvh.sp_hidden);
				};
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
				for (auto sph: _sp_bias_hidden) {
					// h
					inc += hup->get_prob(sph);
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





