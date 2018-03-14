#include "basis_func.hpp"
#include "general.hpp"
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

	IxnParamTraj::IxnParamTraj(std::string name, IxnParamType type, Species *sp, double min, double max, int n, double val0, int n_t) : IxnParamTraj(name,type,sp,nullptr,min,max,n,val0,n_t) {};
	IxnParamTraj::IxnParamTraj(std::string name, IxnParamType type, Species *sp1, Species *sp2, double min, double max, int n, double val0, int n_t) : Grid(name,min,max,n)
	{
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
	IxnParamTraj& IxnParamTraj::operator=(const IxnParamTraj& other) {
		if (this != &other)
		{
			Grid::operator=(other);

			_clean_up();
			_copy(other);
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
	void IxnParamTraj::_clean_up() {
		safeDelArr(_vals);
		safeDelArr(_asleep);
		safeDelArr(_awake);
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
			if (_sp1) {
				if (moment_type==AWAKE) {
					_awake[it] += 1. * _sp1->count() / batch_size;
				} else if (moment_type==ASLEEP) {
					_asleep[it] += 1. * _sp1->count() / batch_size;
				};
			};
		} else if (_type == Jp) {
			if (_sp1 && _sp2) {
				if (moment_type==AWAKE) {
					_awake[it] += 1. * _sp1->nn_count(_sp2) / batch_size;
				} else if (moment_type==ASLEEP) {
					_asleep[it] += 1. * _sp1->nn_count(_sp2) / batch_size;
				};
			};
		};	
	};

	double IxnParamTraj::moments_diff_at_time(int it) {
		return (_awake[it] - _asleep[it]);
	};

	/********************
	Write into an ofstream
	********************/

	void IxnParamTraj::write_vals(std::string dir, int idx, int n_t_traj) const {
		std::ofstream f;
		f.open(dir+name()+"_"+pad_str(idx,4)+".txt");
		for (int i=0; i<n_t_traj; i++) {
			f << _vals[i];
			if (i != n_t_traj-1) { f << "\n"; };
		};
		f.close();
	};
	void IxnParamTraj::write_vals(std::string dir, int idx1, int idx2, int n_t_traj) const {
		std::ofstream f;
		f.open(dir+name()+"_"+pad_str(idx1,4)+"_"+pad_str(idx2,2)+".txt");
		for (int i=0; i<n_t_traj; i++) {
			f << _vals[i];
			if (i != n_t_traj-1) { f << "\n"; };
		};
		f.close();
	};

	void IxnParamTraj::write_moments(std::string dir, int idx, int n_t_traj) const {
		std::ofstream f;
		f.open(dir+name()+"_"+pad_str(idx,4)+".txt");
		for (int i=0; i<n_t_traj; i++) {
			f << _awake[i] << " " << _asleep[i];
			if (i != n_t_traj-1) { f << "\n"; };
		};
		f.close();
	};
	void IxnParamTraj::write_moments(std::string dir, int idx1, int idx2, int n_t_traj) const {
		std::ofstream f;
		f.open(dir+name()+"_"+pad_str(idx1,4)+"_"+pad_str(idx2,2)+".txt");
		for (int i=0; i<n_t_traj; i++) {
			f << _awake[i] << " " << _asleep[i];
			if (i != n_t_traj-1) { f << "\n"; };
		};
		f.close();
	};


};





