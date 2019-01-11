#include "ixn_param_traj.hpp"

// Other headers
#include "../include/dynamicboltz_bits/general.hpp"
#include "species.hpp"
#include "lattice.hpp"
#include "hidden_unit.hpp"
#include "basis_func.hpp"

#include <iostream>
#include <fstream>

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	Ixn Param
	****************************************/

	/********************
	Constructor
	********************/

	IxnParamTraj::IxnParamTraj(std::string name, IxnParamType type, double domain_min, double domain_max, int domain_no_pts, double val0, int n_t, int *t_opt_ptr)
	{
		_name = name;

		_type = type;

		// Domain
		_domain_min = domain_min;
		_domain_max = domain_max;
		_domain_no_pts = domain_no_pts;
		_domain_delta = (_domain_max-_domain_min) / (_domain_no_pts-1);

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

		_t_opt_ptr = t_opt_ptr;

		_is_awake_fixed = false;
		_awake_fixed = nullptr;
	};

	IxnParamTraj::IxnParamTraj(const IxnParamTraj& other) {
		_copy(other);
	};
	IxnParamTraj::IxnParamTraj(IxnParamTraj&& other) {
		_copy(other);
		other._reset();
	};
	IxnParamTraj& IxnParamTraj::operator=(const IxnParamTraj& other) {
		if (this != &other)
		{
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	IxnParamTraj& IxnParamTraj::operator=(IxnParamTraj&& other) {
		if (this != &other)
		{
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
		_name = other._name;

		_type = other._type;

		_domain_min = other._domain_min;
		_domain_max = other._domain_max;
		_domain_delta = other._domain_delta;
		_domain_no_pts = other._domain_no_pts;

		_sp_bias_visible = other._sp_bias_visible;
		_sp_bias_hidden = other._sp_bias_hidden;
		_sp_doublet = other._sp_doublet;
		_sp_triplet = other._sp_triplet;
		_sp_conn = other._sp_conn;

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
		_t_opt_ptr = other._t_opt_ptr;
		if (other._awake_fixed) {
			_awake_fixed = new double[_n_t];
			std::copy( other._awake_fixed, other._awake_fixed + _n_t, _awake_fixed );
		} else {
			_awake_fixed = nullptr;
		};
		_is_awake_fixed = other._is_awake_fixed;
	};
	void IxnParamTraj::_reset()
	{
		_name = "";

		_domain_min = 0.0;
		_domain_max = 0.0;
		_domain_delta = 0.0;
		_domain_no_pts = 0;

		_sp_bias_visible.clear();
		_sp_bias_hidden.clear();
		_sp_doublet.clear();
		_sp_triplet.clear();
		_sp_conn.clear();

		_conns.clear();
		_hidden_units.clear();

		_n_t = 0;
		safeDelArr(_vals);
		_val0 = 0.;
		safeDelArr(_asleep);
		safeDelArr(_awake);
		_bf = nullptr;
		_t_opt_ptr = nullptr;
		if (_awake_fixed) {
			safeDelArr(_awake_fixed);
		};
		_is_awake_fixed = false;
	};
	void IxnParamTraj::_clean_up() {
		safeDelArr(_vals);
		safeDelArr(_asleep);
		safeDelArr(_awake);
		if (_awake_fixed) {
			safeDelArr(_awake_fixed);
		};
	};

	/********************
	Name
	********************/

	std::string IxnParamTraj::get_name() const {
		return _name;
	};

	/********************
	Domain
	********************/

	// Getters
	int IxnParamTraj::get_domain_no_pts() const {
		return _domain_no_pts;
	};
	double IxnParamTraj::get_domain_min() const {
		return _domain_min;
	};
	double IxnParamTraj::get_domain_max() const {
		return _domain_max;
	};
	double IxnParamTraj::get_domain_delta() const {
		return _domain_delta;
	};

	// Get pt in domain
	double IxnParamTraj::get_domain_pt_by_idx(int i) const {
		if (i >= _domain_no_pts) {
			std::cerr << ">>> Error: IxnParamTraj::get_domain_pt_by_idx <<< Idx: " << i << " is out of domain " << _name << " of size " << _domain_no_pts << std::endl;
		};
		return _domain_min + i * _domain_delta;
	};

	// Check if point is in domain
	bool IxnParamTraj::check_if_pt_is_inside_domain(double x) const {
		if (x < _domain_min || x > _domain_max) { 
			return false; 
		} else {
			return true;
		};
	};

	// Get indexes surrounding a point
	// ie point is between i and i+1 where i is returned
	int IxnParamTraj::get_domain_idxs_surrounding_pt(double x) const {
		int i = (x - _domain_min) / _domain_delta;
		if (i==_domain_no_pts-1) {i--;};
		return i;
	};

	// Get fraction of a point between successive points
	double IxnParamTraj::get_domain_frac_between(double x) const {
		return get_domain_frac_between(x,get_domain_idxs_surrounding_pt(x));
	};
	// Second optional specification: the return of the surrounding idxs
	double IxnParamTraj::get_domain_frac_between(double x, int i) const {
		return (x - get_domain_pt_by_idx(i)) / _domain_delta;
	};

	// Print domain range
	void IxnParamTraj::print_domain_range() const {
		std::cout << "Domain: " << _name << " min: " << _domain_min << " max: " << _domain_max << std::endl;
	};

	/********************
	Set time
	********************/

	void IxnParamTraj::set_n_t(int n_t) {
		// Clean
		_clean_up();

		// Set new
		_n_t = n_t;

		// Init new
		_vals = new double[_n_t];
		std::fill_n(_vals,_n_t,0.0);
		_vals[0] = _val0;

		_asleep = new double[_n_t];
		std::fill_n(_asleep,_n_t,0.0);
		_awake = new double[_n_t];
		std::fill_n(_awake,_n_t,0.0);
	};

	/********************
	Add species
	********************/

	void IxnParamTraj::add_species(Species *sp) {
		_sp_bias_visible.push_back(sp);
	};
	void IxnParamTraj::add_species(HiddenSpecies *hsp) {
		_sp_bias_hidden.push_back(hsp);
	};
	void IxnParamTraj::add_species(Species *sp1, Species *sp2) {
		_sp_doublet.push_back(Species2(sp1,sp2));
	};
	void IxnParamTraj::add_species(Species *sp1, Species *sp2, Species *sp3) {
		_sp_triplet.push_back(Species3(sp1,sp2,sp3));
	};
	void IxnParamTraj::add_species(Species *sp, HiddenSpecies *hsp) {
		_sp_conn.push_back(SpeciesVH(sp,hsp));
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
	Set fixed awake
	********************/

	void IxnParamTraj::set_fixed_awake_moment(std::vector<double> vals) {
		if (vals.size() != _n_t) {
			std::cerr << "ERROR - incorrect length for fixed awake moments - needs to be " << _n_t << std::endl;
			exit(EXIT_FAILURE);
		};

		if (!_awake_fixed) {
			_awake_fixed = new double[_n_t];
		};

		for (int i=0; i<vals.size(); i++) {
			_awake_fixed[i] = vals[i];
		};

		_is_awake_fixed = true;
	};

	/********************
	Validation
	********************/

	void IxnParamTraj::validate_setup() const {
		std::cout << "--- Validate ixn param: " << _name << " ---" << std::endl; 
		if (_bf) {
			std::cout << "   Has basis func: " << _bf->name() << std::endl;
		} else {
			std::cerr << "ERROR: no basis func" << std::endl;
			exit(EXIT_FAILURE);
		};
	};

	/********************
	Getters
	********************/

	// At the current time in the optimization
	double IxnParamTraj::get() const {
		return _vals[*_t_opt_ptr];
	};

	// At a specific time
	double IxnParamTraj::get_at_time(int it) const
	{
		return _vals[it];
	};

	/********************
	Get species involved
	********************/

	std::vector<Species*> IxnParamTraj::get_species_bias() const {
		if (_type != Hp) {
			std::cerr << "Error! Not a visible bias." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _sp_bias_visible;
	};
	std::vector<HiddenSpecies*> IxnParamTraj::get_species_hidden_bias() const {
		if (_type != Bp) {
			std::cerr << "Error! Not a hidden bias." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _sp_bias_hidden;
	};
	std::vector<Species2> IxnParamTraj::get_species_doublet() const {
		if (_type != Jp) {
			std::cerr << "Error! Not a J." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _sp_doublet;
	};
	std::vector<Species3> IxnParamTraj::get_species_triplet() const {
		if (_type != Kp) {
			std::cerr << "Error! Not a K." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _sp_triplet;
	};
	std::vector<SpeciesVH> IxnParamTraj::get_species_conn() const {
		if (_type != Wp) {
			std::cerr << "Error! Not a conn W." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _sp_conn;
	};

	/********************
	Calculate the next step
	********************/

	bool IxnParamTraj::calculate_at_time(int it_next, double dt)
	{
		_vals[it_next] = _vals[it_next-1] + dt*_bf->get_at_time(it_next-1);
		return check_if_pt_is_inside_domain(_vals[it_next]);
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
		if (moment_type==AWAKE && _is_awake_fixed) {
			_awake[it] = _awake_fixed[it];
			return;
		};

		if (_type == Hp) {
			for (auto sp: _sp_bias_visible) {
				if (moment_type==AWAKE) {
					_awake[it] += 1. * sp->count() / batch_size;
				} else if (moment_type==ASLEEP) {
					_asleep[it] += 1. * sp->count() / batch_size;
				};
			};
		} else if (_type == Jp) {
			for (auto spd: _sp_doublet) {
				if (moment_type==AWAKE) {
					_awake[it] += 1. * spd.sp1->nn_count(spd.sp2) / batch_size;
				} else if (moment_type==ASLEEP) {
					_asleep[it] += 1. * spd.sp1->nn_count(spd.sp2) / batch_size;
				};
			};
		} else if (_type == Kp) {
			for (auto spt: _sp_triplet) {
				if (moment_type==AWAKE) {
					_awake[it] += 1. * spt.sp1->triplet_count(spt.sp2,spt.sp3) / batch_size;
				} else if (moment_type==ASLEEP) {
					_asleep[it] += 1. * spt.sp1->triplet_count(spt.sp2,spt.sp3) / batch_size;
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
				_awake[it] += 1. * inc / batch_size;
			} else if (moment_type==ASLEEP) {
				_asleep[it] += 1. * inc / batch_size;
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
		std::string fname = dir+_name+"_"+pad_str(idx_opt_step,4);
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
		std::string fname = dir+_name+"_"+pad_str(idx_opt_step,4);
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





