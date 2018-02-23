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

	OptProblem::OptProblem(double nu_min, double nu_max, int n_nu, double t_max, int n_t, double nu_init, int batch_size, int n_annealing, int box_length, double dopt, int n_opt) : _f(nu_min,nu_max,n_nu), _latt(box_length)
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
		this->_dopt = dopt;
		this->_n_opt = n_opt;

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
	};

	/********************
	Set properties
	********************/

	void OptProblem::add_species(std::string sp) {
		_species.push_back(sp);

		// Add to the lattice
		_latt.add_species(&(_species.back()));

		// Add entry for the moments
		_moms_awake[&(_species.back())] = 0.0;
		_moms_asleep[&(_species.back())] = 0.0;
	};

	void OptProblem::add_fname(std::string f) {
		_fnames.push_back(f);
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

	void OptProblem::solve(bool verbose) {
		
		// Write the moments
		std::ofstream fmoments;

		// Filename to read
		std::string fname;

		// Map of interactions
		std::map<std::string,double> h_dict;
		std::map<std::string,std::map<std::string,double>> j_dict;

		// BFs update
		double *df = new double[_n_nu];

		// Iterate over optimization steps
		for (int i_opt=0; i_opt<_n_opt; i_opt++)
		{
			std::cout << "Opt step " << i_opt << " / " << _n_opt-1 << std::endl;

			// Write the current basis funcs
			write_bfs("data/F/"+pad_str(i_opt,4)+".txt");

			// Solve the current nu trajectory
			solve_nu_traj();
			if (verbose) {
				for (int i=0; i<_n_t; i++) {
					std::cout << _nu_traj[i] << " ";
				};
				std::cout << std::endl;
			};

			// Write the nu trajectory
			write_nu_traj("data/nu_traj/"+pad_str(i_opt,4)+".txt");

			// Solve the variational problem traj
			solve_var_traj();

			// Write the nu trajectory
			write_var_traj("data/var_traj/"+pad_str(i_opt,4)+".txt");

			// Pick a random batch
			std::vector<std::string> fnames, fnames_possible=_fnames;
			std::vector<std::string>::iterator itf;
			for (int i_batch=0; i_batch<_n_batch; i_batch++) {
				itf = fnames_possible.begin();
				std::advance(itf,randI(0,fnames_possible.size()-1));
				fnames.push_back(*itf);
				fnames_possible.erase(itf);
			};

			// Reset update
			std::fill_n(df, _n_nu, 0.);

			// Get ready to write the moments
			fmoments.open("data/moments/"+pad_str(i_opt,4)+".txt");

			// Go through all times
			for (int i_time=0; i_time < _n_t_soln; i_time++)
			{
				if (verbose) {
					std::cout << "time: " << i_time << std::flush;
				};

				// Reset the moments
				for (auto it=_species.begin(); it!=_species.end(); it++)
				{
					_moms_awake[&(*it)] = 0.0;
					_moms_asleep[&(*it)] = 0.0;
				};

				// Write the interactions at this timestep into the species
				for (auto it1=_species.begin(); it1!=_species.end(); it1++)
				{
					it1->h = _nu_traj[i_time];
					for (auto it2=_species.begin(); it2!=_species.end(); it2++)
					{
						it1->j[&(*it2)] = 0.0;
					};
				};

				// Go through all samples in the batch
				for (int i_batch=0; i_batch<_n_batch; i_batch++) 
				{
					if (verbose) {
						std::cout << "." << std::flush;
					};

					// Read in batch at this timestep
					fname = fnames[i_batch] + pad_str(i_time,4) + ".txt";
					_latt.read_from_file(fname);

					// Record the moments
					for (auto it=_species.begin(); it!=_species.end(); it++)
					{
						_moms_awake[&(*it)] += 1.0 * it->count / _n_batch;
					};

					// Anneal
					_latt.anneal(_n_annealing);

					// Record the moments
					for (auto it=_species.begin(); it!=_species.end(); it++)
					{
						_moms_asleep[&(*it)] += 1.0 * it->count / _n_batch;
					};
				};

				if (verbose) {
					std::cout << std::endl << _moms_awake[&(_species.front())] << " " << _moms_asleep[&(_species.front())] << std::endl;
				};

				// Write the moments
				fmoments << _t_grid[i_time] << " " << _moms_awake[&(_species.front())] << " " << _moms_asleep[&(_species.front())] << "\n";

				// Store the update AND DAMPEN
				for (int j=0; j<_n_nu; j++) {
					df[j] -= _dopt * _dt * (_moms_asleep[&(_species.front())] - _moms_awake[&(_species.front())]) * _var_traj[i_time][j];//* exp(-pow(_nu_traj[i_time]-_nu_grid[j],2)/(2.*0.1));
				};
			};

			// Close the moments file
			fmoments.close();

			// Update the basis functions
			_f.update(df);
		};

		// Write the final basis funcs once more
		write_bfs("data/F/"+pad_str(_n_opt,4)+".txt");

		// Free update
		delete[] df;
	};

	/****************************************
	OptProblem - PRIVATE
	****************************************/

	/********************
	Delta function at some time index, nu value
	********************/

	double OptProblem::_delta(double nu1, double nu2) 
	{
		return exp(-pow(nu1-nu2,2)/(2.*0.1)) / sqrt(2. * M_PI * 0.1);
	};

};

