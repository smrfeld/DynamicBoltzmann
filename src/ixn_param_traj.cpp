#include "basis_func.hpp"
#include "../include/general.hpp"
#include <iostream>
#include <fstream>

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

	IxnParamTraj::IxnParamTraj(std::string name, IxnParamType type, Species *sp, double min, double max, int n, double val0, int n_t) : IxnParamTraj(name,type,sp,nullptr,min,max,n,val0,n_t) {
		if (type != Wp && type != Hp && type != Bp) {
			std::cerr << "ERROR: wrong IxnParamTraj constructor?" << std::endl;
			exit(EXIT_FAILURE);
		};
	};
	IxnParamTraj::IxnParamTraj(std::string name, IxnParamType type, Species *sp1, Species *sp2, double min, double max, int n, double val0, int n_t) : Grid(name,min,max,n)
	{
		if (type == Jp && (sp1 == nullptr || sp2==nullptr)) {
			std::cerr << "ERROR: wrong IxnParamTraj constructor?" << std::endl;
			exit(EXIT_FAILURE);
		};

		_type = type;
		_sp1 = sp1;
		_sp2 = sp2;
		_val0 = val0;
		_n_t = n_t;

		_vals = new double[_n_t];
		std::fill_n(_vals,_n_t,0.0);
		_vals[0] = _val0;

		_asleep = new double[_n_t];
		std::fill_n(_asleep,_n_t,0.0);
		_awake = new double[_n_t];
		std::fill_n(_awake,_n_t,0.0);

		_bf = nullptr;
	};

	IxnParamTraj::IxnParamTraj(const IxnParamTraj& other) : Grid(other) {
		_copy(other);
	};
	IxnParamTraj::IxnParamTraj(IxnParamTraj&& other) : Grid(other) {
		_copy(other);
		other._reset();
	};
	IxnParamTraj& IxnParamTraj::operator=(const IxnParamTraj& other) {
		if (this != &other)
		{
			Grid::operator=(other);

			_clean_up();
			_copy(other);
		};
		return *this;
	};
	IxnParamTraj& IxnParamTraj::operator=(IxnParamTraj&& other) {
		if (this != &other)
		{
			Grid::operator=(other);

			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
	};
	IxnParamTraj::~IxnParamTraj() {
		_clean_up();
	};
	void IxnParamTraj::_copy(const IxnParamTraj& other)
	{
		_type = other._type;
		_sp1 = other._sp1;
		_sp2 = other._sp2;
		_conns = other._conns;
		_hidden_units = other._hidden_units;
		_n_t = other._n_t;
		_vals = new double[_n_t];
		std::copy( other._vals, other._vals + _n_t, _vals );
		_val0 = other._val0;
		_asleep = new double[_n_t];
		std::copy( other._asleep, other._asleep + _n_t, _asleep );
		_awake = new double[_n_t];
		std::copy( other._awake, other._awake + _n_t, _awake );
		_bf = other._bf;
	};
	void IxnParamTraj::_reset()
	{
		_sp1 = nullptr;
		_sp2 = nullptr;
		_conns.clear();
		_hidden_units.clear();
		_n_t = 0;
		safeDelArr(_vals);
		_val0 = 0.;
		safeDelArr(_asleep);
		safeDelArr(_awake);
		_bf = nullptr;
	};
	void IxnParamTraj::_clean_up() {
		safeDelArr(_vals);
		safeDelArr(_asleep);
		safeDelArr(_awake);
	};

	/********************
	Check if this ixn param is...
	********************/

	bool IxnParamTraj::is_h_with_species(std::string species_name) const {
		if (_type == Hp && _sp1->name() == species_name) {
			return true;
		} else {
			return false;
		};
	};

	bool IxnParamTraj::is_w_with_species(std::string species_name) const {
		if (_type == Wp && _sp1->name() == species_name) {
			return true;
		} else {
			return false;
		};
	};

	bool IxnParamTraj::is_j_with_species(std::string species_name_1, std::string species_name_2) const {
		if (_type == Jp) {
			if ((_sp1->name() == species_name_1 && _sp2->name() == species_name_2) || (_sp1->name() == species_name_2 && _sp2->name() == species_name_1)) {
				return true;
			};
		};
		return false;
	};

	bool IxnParamTraj::is_b_with_species(std::string species_name) const {
		if (_type == Bp && _sp1->name() == species_name) {
			return true;
		} else {
			return false;
		};
	};

	/********************
	If Wp;
	Add a visible->hidden unit connection
	********************/

	void IxnParamTraj::add_visible_hidden_connection(Site *sptr, HiddenUnit *hup) {
		if (_type != Wp) {
			std::cerr << "ERROR: adding visible to hidden connection to wrong type of IxnParam" << std::endl;
			exit(EXIT_FAILURE);
		};
		_conns.push_back(std::make_pair(sptr,hup));
	};

	/********************
	If Bp;
	Add a hidden unit to monitor
	********************/

	void IxnParamTraj::add_hidden_unit(HiddenUnit *hup) {
		if (_type != Bp) {
			std::cerr << "ERROR: adding visible to hidden connection to wrong type of IxnParam" << std::endl;
			exit(EXIT_FAILURE);
		};
		_hidden_units.push_back(hup);
	};

	/********************
	Set/check basis func pointer
	********************/

	void IxnParamTraj::set_basis_func_ptr(BasisFunc* bf) {
		_bf = bf;
	};
	bool IxnParamTraj::is_bf(BasisFunc *bf) {
		if (_bf==bf) { 
			return true;
		} else { 
			return false; 
		};
	};

	/********************
	Set IC
	********************/

	void IxnParamTraj::set_init_cond(double val) {
		_val0 = val;
		_vals[0] = _val0;
	};

	/********************
	Validation
	********************/

	void IxnParamTraj::validate_setup() const {
		std::cout << "--- Validate ixn param: " << name() << " ---" << std::endl; 
		if (_bf) {
			std::cout << "   Has basis func: " << _bf->name() << std::endl;
		} else {
			std::cerr << "ERROR: no basis func" << std::endl;
			exit(EXIT_FAILURE);
		};
	};

	/********************
	Getters/setters
	********************/

	double IxnParamTraj::get_at_time(int it) const
	{
		return _vals[it];
	};

	/********************
	Calculate the next step
	********************/

	bool IxnParamTraj::calculate_at_time(int it_next, double dt)
	{
		_vals[it_next] = _vals[it_next-1] + dt*_bf->get_at_time(it_next-1);
		return in_grid(_vals[it_next]);
	};

	/********************
	Moments from lattice
	********************/

	void IxnParamTraj::moments_reset() 
	{
		for (int it=0; it<_n_t; it++) {
			_asleep[it] = 0.;
			_awake[it] = 0.;
		};
	};
	void IxnParamTraj::moments_retrieve_at_time(MomentType moment_type, int it) {
		moments_retrieve_at_time(moment_type, it, 1);
	};
	void IxnParamTraj::moments_retrieve_at_time(MomentType moment_type, int it, int batch_size)
	{
		if (_type == Hp) {
			if (moment_type==AWAKE) {
				_awake[it] += 1. * _sp1->count() / batch_size;
			} else if (moment_type==ASLEEP) {
				_asleep[it] += 1. * _sp1->count() / batch_size;
			};
		} else if (_type == Jp) {
			if (moment_type==AWAKE) {
				_awake[it] += 1. * _sp1->nn_count(_sp2) / batch_size;
			} else if (moment_type==ASLEEP) {
				_asleep[it] += 1. * _sp1->nn_count(_sp2) / batch_size;
			};
		} else if (_type == Wp) {
			double inc = 0.0;
			// Run through all connections
			for (auto c: _conns) {
				if (c.first->sp == _sp1) { // site occupied with the correct species
					inc += c.second->get(); // v * h
				};
			};
			if (moment_type==AWAKE) {
				_awake[it] += inc / batch_size;
			} else if (moment_type==ASLEEP) {
				_asleep[it] += inc / batch_size;
			};
		} else if (_type == Bp) {
			// Run through all hidden units that are monitored
			double inc = 0.0;
			for (auto hup: _hidden_units) {
				inc += hup->get();
			};
			if (moment_type==AWAKE) {
				_awake[it] += 1. * inc / batch_size;
			} else if (moment_type==ASLEEP) {
				_asleep[it] += 1. * inc / batch_size;
			};
		};
	};

	double IxnParamTraj::moments_diff_at_time(int it) {
		return (_awake[it] - _asleep[it]);
	};
	double IxnParamTraj::moments_get_at_time(MomentType moment_type, int it) const {
		if (moment_type == AWAKE) {
			return _awake[it];
		} else {
			return _asleep[it];
		};
	};

	/********************
	Write into an ofstream
	********************/

	void IxnParamTraj::write_vals(std::string dir, int idx_opt_step, std::vector<int> idxs, int n_t_traj) const {
		std::ofstream f;
		std::string fname = dir+name()+"_"+pad_str(idx_opt_step,4);
		for (auto idx: idxs) {
			fname += "_" + pad_str(idx,4);
		};
		f.open(fname+".txt");
		for (int i=0; i<n_t_traj; i++) {
			f << _vals[i];
			if (i != n_t_traj-1) { f << "\n"; };
		};
		f.close();
	};

	void IxnParamTraj::write_moments(std::string dir, int idx_opt_step, std::vector<int> idxs, int n_t_traj) const {
		std::ofstream f;
		std::string fname = dir+name()+"_"+pad_str(idx_opt_step,4);
		for (auto idx: idxs) {
			fname += "_" + pad_str(idx,4);
		};
		f.open(fname+".txt");
		for (int i=0; i<n_t_traj; i++) {
			f << _awake[i] << " " << _asleep[i];
			if (i != n_t_traj-1) { f << "\n"; };
		};
		f.close();
	};


};





