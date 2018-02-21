#include "dynamic_boltzmann.hpp"
#include <numeric> 
#include "general.hpp"
#include "math.h"
#include <fstream>

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	OptProblem
	****************************************/

	/********************
	Constructor
	********************/

	OptProblem::OptProblem(double nu_min, double nu_max, int n_nu, double t_max, int n_t, double nu_init, int batch_size, int n_annealing, std::vector<std::string> species_list, int box_length, double dopt) : _f(nu_min,nu_max,n_nu)
	{
		this->_t_max = t_max;
		this->_n_t = n_t;
		this->_n_t_soln = n_t;
		this->_dt = t_max / (n_t - 1);
		this->_nu_init = nu_init;
		this->_nu_min = nu_min;
		this->_nu_max = nu_max;
		this->_n_nu = n_nu;
		this->_dnu = (nu_max-nu_min) / (n_nu - 1);
		this->_n_batch = batch_size;
		this->_n_annealing = n_annealing;
		this->_box_length = box_length;
		this->_species = species_list;
		this->_dopt = dopt;

		// Species
		for (auto s: species_list) {
			this->_moms_awake[s] = 0.0;
			this->_moms_asleep[s] = 0.0;
		};

		// Data
		this->_data = new Lattice[batch_size];
		for (int i=0; i<batch_size; i++) {
			this->_data[i] = Lattice(box_length);
		};

		// Fill the _nu_grid, _t_grid incrementally	
		this->_nu_grid = new double[n_nu];
	    for (int i=0; i<n_nu; i++) {
	    	this->_nu_grid[i] = nu_min + i*this->_dnu;
	    };
		this->_t_grid = new double[n_t];
	    for (int i=0; i<n_t; i++) {
	    	this->_t_grid[i] = i*this->_dt;
	    };

	    // Initialize nu traj, var traj
	    this->_nu_traj = new double[n_t];
		std::fill_n(this->_nu_traj, n_t, 0.0);
	    this->_var_traj = new double*[n_t];
	    for (int i=0; i<n_t; i++) {
	    	this->_var_traj[i] = new double[n_nu];
		};
	};
	OptProblem::~OptProblem() {
		safeDelArr(this->_nu_traj);
		safeDelArr(this->_t_grid);
		safeDelArr(this->_nu_grid);
		if (this->_var_traj != nullptr)
		{
			for (int i=0; i<_n_t; i++) {
				safeDelArr(this->_var_traj[i]);
			};
			delete[] this->_var_traj;
		};
		safeDelArr(this->_data);
	};

	/********************
	Solve for nu traj
	********************/

	void OptProblem::solve_nu_traj() 
	{
		_nu_traj[0] = _nu_init;
		_t_grid[0] = 0.;
		double n;
		for (int i=0; i<_n_t-1; i++) {
			n = _nu_traj[i] + _dt * _f.get_val(_nu_traj[i]);
			if (n <= _nu_max && n >= _nu_min) {
				_nu_traj[i+1] = n;
				_t_grid[i+1] = (i+1)*_dt;
			} else {
				_n_t_soln = i+1;
				return;
			};
		};
		_n_t_soln = _n_t;
	};

	/********************
	Solve for var traj
	********************/

	void OptProblem::solve_var_traj()
	{
		// Initial
		for (int j=0; j<_n_nu; j++) {
			_var_traj[0][j] = 0.;
		};

		// Go through all times
		for (int i=0; i<_n_t_soln-1; i++) {
			// Go through all nu
			for (int j=0; j<_n_nu; j++) {
				_var_traj[i+1][j] = _var_traj[i][j] + _dt * ( _delta(_nu_traj[i],_nu_grid[j]) + _f.get_deriv(_nu_traj[i]) * _var_traj[i][j] );
			};
		};
	};

	/********************
	Print
	********************/

	void OptProblem::print_nu_traj() {
		for (int i=0; i<_n_t_soln; i++) 
		{ 
			std::cout << _t_grid[i] << " " << _nu_traj[i] << std::endl;
		};
	};
	void OptProblem::print_var_traj() {
		for (int i=0; i<_n_t_soln; i++) 
		{
			for (int j=0; j<_n_nu; j++) {
				std::cout << _t_grid[i] << " " << _nu_grid[j] << " " << _var_traj[i][j] << std::endl;
			};
		};
	};

	/********************
	Write
	********************/

	void OptProblem::write_nu_traj(std::string fname) {
		std::ofstream f;
		f.open (fname);
		for (int i=0; i<_n_t_soln; i++) {
			f << _t_grid[i] << " " << _nu_traj[i];
			if (i != _n_t_soln-1) { f << "\n";};
		};
		f.close();
	};
	void OptProblem::write_var_traj(std::string fname) {
		std::ofstream f;
		f.open (fname);
		for (int i=0; i<_n_t_soln; i++) {
			for (int j=0; j<_n_nu; j++) {
				f << _t_grid[i] << " " << _nu_grid[j] << " " << _var_traj[i][j];
				if (i != _n_t_soln-1 || j != _n_nu-1) { f << "\n";};
			};
		};
		f.close();
	};
	void OptProblem::write_bfs(std::string fname) {
		_f.write_to_file(fname);
	};

	/********************
	Run main optimization loop
	********************/

	void OptProblem::solve() {

		// Solve the current nu trajectory
		solve_nu_traj();

		// Solve the variational problem traj
		solve_var_traj();

		// Pick random batch indexes
		// ...

		// Map of interactions
		std::map<std::string,double> h_dict;
		std::map<StrPair,double> j_dict;

		// BFs update
		double *df = new double[_n_nu];
		std::fill_n(df, _n_nu, 0.);

		// Go through all times
		for (int i_time=0; i_time < _n_t_soln; i_time++)
		{
			std::cout << "time: " << i_time << std::endl;

			// Read in batch at this timestep
			for (int i_batch=0; i_batch<_n_batch; i_batch++) {
				_data[i_batch].read_from_file("latt_st_test.txt");

				// Record the moments
				for (auto sp: _species) {
					_moms_awake[sp] += _data[i_batch].get_count(sp) / _n_batch;
				};
			};

			// Form the dictionary of interactions at this timestep
			for (int is1=0; is1<_species.size(); is1++) {
				h_dict[_species[is1]] = _nu_traj[i_time];
				for (int is2=is1; is2<_species.size(); is2++) {
					j_dict[StrPair(_species[is1],_species[is2])] = 0.0;
				};
			};

			// Anneal the batch
			for (int i_batch=0; i_batch<_n_batch; i_batch++) {
				_data[i_batch].anneal(h_dict,j_dict,_n_annealing);

				// Record the moments
				for (auto sp: _species) {
					_moms_asleep[sp] += _data[i_batch].get_count(sp) / _n_batch;
				};
			};

			// Store the update
			for (int j=0; j<_n_nu; j++) {
				df[j] -= _dopt * _dt * (_moms_asleep[_species[0]] - _moms_awake[_species[0]]) * _var_traj[i_time][j];
			};
		};

		// Update the basis functions
		_f.update(df);
	};

	/****************************************
	OptProblem - PRIVATE
	****************************************/

	/********************
	Delta function at some time index, nu value
	********************/

	double OptProblem::_delta(double nu1, double nu2) 
	{
		return exp(-pow(nu1-nu2,2)/(2.*0.01)) / sqrt(2. * M_PI * 0.01);
	};

};

