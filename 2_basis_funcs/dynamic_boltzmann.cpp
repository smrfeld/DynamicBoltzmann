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

	OptProblem::OptProblem(double h_min, double h_max, int n_h, double j_min, double j_max, int n_j, double t_max, int n_t, double h_init, double j_init, int batch_size, int n_annealing, int box_length, double dopt, int n_opt) : _f_h(h_min,h_max,n_h,j_min,j_max,n_j), _f_j(h_min,h_max,n_h,j_min,j_max,n_j), _latt(box_length)
	{
		// Times
		this->_t_max = t_max;
		this->_n_t = n_t;
		this->_n_t_soln = 0;
		this->_dt = t_max / (n_t - 1);

		// h, j
		this->_h_init = h_init;
		this->_j_init = j_init;
		this->_h_min = h_min;
		this->_h_max = h_max;
		this->_n_h = n_h;
		this->_dh = (h_max-h_min) / (n_h - 1);
		this->_j_min = j_min;
		this->_j_max = j_max;
		this->_n_j = n_j;
		this->_dj = (j_max-j_min) / (n_j - 1);

		// Optimization params
		this->_n_batch = batch_size;
		this->_n_annealing = n_annealing;
		this->_box_length = box_length;
		this->_dopt = dopt;
		this->_n_opt = n_opt;

		// Fill the _h_grid, _j_grid, _t_grid incrementally	
		this->_h_grid = new double[n_h];
	    for (int i=0; i<n_h; i++) {
	    	this->_h_grid[i] = h_min + i*_dh;
	    };
		this->_j_grid = new double[n_j];
	    for (int i=0; i<n_j; i++) {
	    	this->_j_grid[i] = j_min + i*_dj;
	    };
		this->_t_grid = new double[n_t];
	    for (int i=0; i<n_t; i++) {
	    	this->_t_grid[i] = i*_dt;
	    };

	    // Initialize nu traj, var traj
	    this->_h_traj = new double[n_t];
		std::fill_n(this->_h_traj, n_t, 0.0);
	    this->_j_traj = new double[n_t];
		std::fill_n(this->_j_traj, n_t, 0.0);

	    this->_var_hh_traj = new double**[n_t];
	    this->_var_hj_traj = new double**[n_t];
	    this->_var_jh_traj = new double**[n_t];
	    this->_var_jj_traj = new double**[n_t];
	    for (int i=0; i<n_t; i++) {
		    this->_var_hh_traj[i] = new double*[n_h];
		    this->_var_hj_traj[i] = new double*[n_h];
		    this->_var_jh_traj[i] = new double*[n_h];
		    this->_var_jj_traj[i] = new double*[n_h];
		    for (int j=0; j<n_h; j++) {
			    this->_var_hh_traj[i][j] = new double[n_j];
			    this->_var_hj_traj[i][j] = new double[n_j];
			    this->_var_jh_traj[i][j] = new double[n_j];
			    this->_var_jj_traj[i][j] = new double[n_j];
	    		std::fill_n(this->_var_hh_traj[i][j], n_j, 0.0);
				std::fill_n(this->_var_hj_traj[i][j], n_j, 0.0);
				std::fill_n(this->_var_jh_traj[i][j], n_j, 0.0);
				std::fill_n(this->_var_jj_traj[i][j], n_j, 0.0);
		    };
		};
	};
	OptProblem::~OptProblem() {
		safeDelArr(_h_traj);
		safeDelArr(_j_traj);
		safeDelArr(_t_grid);
		safeDelArr(_h_grid);
		safeDelArr(_j_grid);
		if (_var_hh_traj != nullptr)
		{
			for (int i=0; i<_n_t; i++) {
				for (int j=0; j<_n_h; j++) {
					safeDelArr(_var_hh_traj[i][j]);
				};
				delete[] _var_hh_traj[i];
			};
			delete[] _var_hh_traj;
		};
		if (_var_hj_traj != nullptr)
		{
			for (int i=0; i<_n_t; i++) {
				for (int j=0; j<_n_h; j++) {
					safeDelArr(_var_hj_traj[i][j]);
				};
				delete[] _var_hj_traj[i];
			};
			delete[] _var_hj_traj;
		};
		if (_var_jh_traj != nullptr)
		{
			for (int i=0; i<_n_t; i++) {
				for (int j=0; j<_n_h; j++) {
					safeDelArr(_var_jh_traj[i][j]);
				};
				delete[] _var_jh_traj[i];
			};
			delete[] _var_jh_traj;
		};
		if (_var_jj_traj != nullptr)
		{
			for (int i=0; i<_n_t; i++) {
				for (int j=0; j<_n_h; j++) {
					safeDelArr(_var_jj_traj[i][j]);
				};
				delete[] _var_jj_traj[i];
			};
			delete[] _var_jj_traj;
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
		_moms_h_asleep[&(_species.back())] = 0.0;
		_moms_h_awake[&(_species.back())] = 0.0;
		for (auto it=_species.begin(); it!=_species.end(); it++) {
			_moms_j_asleep[&(_species.back())][&(*it)] = 0.0;
			_moms_j_asleep[&(*it)][&(_species.back())] = 0.0;
			_moms_j_awake[&(_species.back())][&(*it)] = 0.0;
			_moms_j_awake[&(*it)][&(_species.back())] = 0.0;
		};
	};

	void OptProblem::add_fname(std::string f) {
		_fnames.push_back(f);
	};

	/********************
	Solve for nu traj
	********************/

	void OptProblem::solve_hj_traj() 
	{
		_h_traj[0] = _h_init;
		_j_traj[0] = _j_init;
		_t_grid[0] = 0.;
		double hnew, jnew;
		for (int i=0; i<_n_t-1; i++) {
			hnew = _h_traj[i] + _dt * _f_h.get_val(_h_traj[i],_j_traj[i]);
			jnew = _j_traj[i] + _dt * _f_j.get_val(_h_traj[i],_j_traj[i]);
			if (hnew <= _h_max && hnew >= _h_min && jnew <= _j_max && jnew >= _j_min) {
				_h_traj[i+1] = hnew;
				_j_traj[i+1] = jnew;
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
		for (int i=0; i<_n_h; i++) {
			for (int j=0; j<_n_j; j++) {
				_var_hh_traj[0][i][j] = 0.;
				_var_hj_traj[0][i][j] = 0.;
				_var_jh_traj[0][i][j] = 0.;
				_var_jj_traj[0][i][j] = 0.;
			};
		};

		// Go through all times
		for (int t=0; t<_n_t_soln-1; t++) {
			// Go through all h,j
			for (int i=0; i<_n_h; i++) {
				for (int j=0; j<_n_j; j++) {
					_var_hh_traj[t+1][i][j] = _var_hh_traj[t+1][i][j] + _dt * (
						_delta(_h_traj[t],_h_grid[i],_j_traj[t],_j_grid[j])
						+
						_f_h.get_deriv(_h_traj[t],_j_traj[t], true, false) * _var_hh_traj[t][i][j]
						+
						_f_h.get_deriv(_h_traj[t],_j_traj[t], false, true) * _var_jh_traj[t][i][j]
						);
					_var_hj_traj[t+1][i][j] = _var_hh_traj[t+1][i][j] + _dt * (
						_f_h.get_deriv(_h_traj[t],_j_traj[t], true, false) * _var_hh_traj[t][i][j]
						+
						_f_h.get_deriv(_h_traj[t],_j_traj[t], false, true) * _var_jj_traj[t][i][j]
						);
					_var_jh_traj[t+1][i][j] = _var_hh_traj[t+1][i][j] + _dt * (
						_delta(_h_traj[t],_h_grid[i],_j_traj[t],_j_grid[j])
						+
						_f_j.get_deriv(_h_traj[t],_j_traj[t], true, false) * _var_hh_traj[t][i][j]
						+
						_f_j.get_deriv(_h_traj[t],_j_traj[t], false, true) * _var_jh_traj[t][i][j]
						);
					_var_jj_traj[t+1][i][j] = _var_hh_traj[t+1][i][j] + _dt * (
						_delta(_h_traj[t],_h_grid[i],_j_traj[t],_j_grid[j])
						+
						_f_j.get_deriv(_h_traj[t],_j_traj[t], true, false) * _var_hj_traj[t][i][j]
						+
						_f_j.get_deriv(_h_traj[t],_j_traj[t], false, true) * _var_jj_traj[t][i][j]
						);
				};
			};
		};
	};

	/********************
	Write
	********************/

	void OptProblem::write_hj_traj(std::string fname) {
		std::ofstream f;
		f.open (fname);
		for (int i=0; i<_n_t_soln; i++) {
			f << _t_grid[i] << " " << _h_traj[i] << " " << _j_traj[i];
			if (i != _n_t_soln-1) { f << "\n";};
		};
		f.close();
	};
	void OptProblem::write_var_traj(std::string fname) {
		std::ofstream f;
		f.open (fname);
		for (int t=0; t<_n_t_soln; t++) {
			for (int i=0; i<_n_h; i++) {
				for (int j=0; j<_n_j; j++) {
					f << _t_grid[t] << " " << _h_grid[i] << " " << _j_grid[j] << " " << _var_hh_traj[t][i][j] << " " << _var_hj_traj[t][i][j] << " " << _var_jh_traj[t][i][j] << " " << _var_jj_traj[t][i][j];
					if (t != _n_t_soln-1 || i != _n_h-1 || j != _n_j-1) { f << "\n";};
				};
			};
		};
		f.close();
	};
	void OptProblem::write_bfs(std::string fname_h, std::string fname_j) {
		_f_h.write_to_file(fname_h);
		_f_j.write_to_file(fname_j);
	};

	/********************
	Run main optimization loop
	********************/

	void OptProblem::solve(bool verbose) {
		
		// Write the moments
		std::ofstream fmoments;

		// Filename to read
		std::string fname;

		// BFs update
		double **df_h = new double*[_n_h];
		double **df_j = new double*[_n_h];
		for (int i=0; i<_n_h; i++) {
			df_h[i] = new double[_n_j];
			df_j[i] = new double[_n_j];
		};

		// Iterate over optimization steps
		for (int i_opt=0; i_opt<_n_opt; i_opt++)
		{
			std::cout << "Opt step " << i_opt << " / " << _n_opt-1 << std::endl;

			// Write the current basis funcs
			write_bfs("data/Fh/"+pad_str(i_opt,4)+".txt","data/Fj/"+pad_str(i_opt,4)+".txt");

			// Solve the current nu trajectory
			solve_hj_traj();
			if (verbose) {
				for (int i=0; i<_n_t; i++) {
					std::cout << _h_traj[i] << " ";
				};
				std::cout << std::endl;
			};

			// Write the nu trajectory
			write_hj_traj("data/traj/"+pad_str(i_opt,4)+".txt");

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
			for (int i=0; i<_n_h; i++) {
				std::fill_n(df_h[i], _n_j, 0.);
				std::fill_n(df_j[i], _n_j, 0.);
			};

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
					_moms_h_awake[&(*it)] = 0.0;
					_moms_h_asleep[&(*it)] = 0.0;
					for (auto jt=_species.begin(); jt!=_species.end(); jt++)
					{
						_moms_j_awake[&(*it)][&(*jt)] = 0.0;
						_moms_j_asleep[&(*it)][&(*jt)] = 0.0;
					};
				};

				// Write the interactions at this timestep into the species
				for (auto it1=_species.begin(); it1!=_species.end(); it1++)
				{
					it1->h = _h_traj[i_time];
					for (auto it2=_species.begin(); it2!=_species.end(); it2++)
					{
						it1->j[&(*it2)] = _j_traj[i_time];
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
						_moms_h_awake[&(*it)] += 1.0 * it->count / _n_batch;
						for (auto jt=it; jt!=_species.end(); jt++)
						{
							_moms_j_awake[&(*it)][&(*jt)] += 1.0 * it->nn_count[&(*jt)] / _n_batch;
						};
					};

					// Anneal
					_latt.anneal(_n_annealing);

					// Record the moments
					for (auto it=_species.begin(); it!=_species.end(); it++)
					{
						_moms_h_asleep[&(*it)] += 1.0 * it->count / _n_batch;
						for (auto jt=it; jt!=_species.end(); jt++)
						{
							_moms_j_asleep[&(*it)][&(*jt)] += 1.0 * it->nn_count[&(*jt)] / _n_batch;
						};
					};
				};

				if (verbose) {
					std::cout << std::endl << _moms_h_awake[&(_species.front())] << " " << _moms_h_asleep[&(_species.front())] << " " << _moms_j_awake[&(_species.front())][&(_species.front())] << " " << _moms_j_asleep[&(_species.front())][&(_species.front())] << std::endl;
				};

				// Write the moments
				fmoments << _t_grid[i_time] << " " << _moms_h_awake[&(_species.front())] << " " << _moms_h_asleep[&(_species.front())] << " " << _moms_j_awake[&(_species.front())][&(_species.front())] << " " << _moms_j_asleep[&(_species.front())][&(_species.front())] << "\n";

				// Store the update
				for (int i=0; i<_n_h; i++) {
					for (int j=0; j<_n_j; j++) {
						df_h[i][j] -= _dopt * _dt * ( 
							(_moms_h_asleep[&(_species.front())] - _moms_h_awake[&(_species.front())]) * _var_hh_traj[i_time][i][j] 
							+
							(_moms_j_asleep[&(_species.front())][&(_species.front())] - _moms_j_awake[&(_species.front())][&(_species.front())]) * _var_jh_traj[i_time][i][j] 
							);
						df_j[i][j] -= _dopt * _dt * ( 
							(_moms_h_asleep[&(_species.front())] - _moms_h_awake[&(_species.front())]) * _var_hj_traj[i_time][i][j] 
							+
							(_moms_j_asleep[&(_species.front())][&(_species.front())] - _moms_j_awake[&(_species.front())][&(_species.front())]) * _var_jj_traj[i_time][i][j] 
							);
					};
				};
			};

			// Close the moments file
			fmoments.close();

			// Update the basis functions
			_f_h.update(df_h);
			_f_j.update(df_j);
		};

		// Write the final basis funcs once more
		write_bfs("data/Fh/"+pad_str(_n_opt,4)+".txt","data/Fj/"+pad_str(_n_opt,4)+".txt");

		// Free update
		for (int i=0; i<_n_h; i++) {
			safeDelArr(df_h[i]);
			safeDelArr(df_j[i]);
		};
		delete[] df_h; 
		delete[] df_j;
	};

	/****************************************
	OptProblem - PRIVATE
	****************************************/

	/********************
	Delta function at some time index, nu value
	********************/

	double OptProblem::_delta(double h1, double h2, double j1, double j2)
	{
		return exp(-pow(h1-h2,2)/(2.*0.1)) * exp(-pow(j1-j2,2)/(2.*0.1)) / (2. * M_PI * 0.1);
	};

};

