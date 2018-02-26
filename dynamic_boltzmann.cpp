#include "dynamic_boltzmann.hpp"
#include <numeric> 
#include "general.hpp"
#include "math.h"
#include <fstream>
#include <iostream>

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

	OptProblem::OptProblem(int n_dim, double t_max, int n_t, int batch_size, int n_annealing, int box_length, double dopt, int n_opt, std::vector<Dim> dims, std::vector<double> init, std::vector<std::string> species) : _latt(box_length), _time("time",0.0,t_max,n_t)
	{
		if (DIAG_INIT) {
			std::cout << "Initializing..." << std::flush;
		};

		// Dimensions
		_n_dim = n_dim;
		_dims = dims;

		// Initial points
		_init = init;

		// Initialize basis funcs
		if (DIAG_INIT) {
			std::cout << "basis funcs..." << std::flush;
		};
		for (int id=0; id<_n_dim; id++) {
			// Basis functions
			_bfs.push_back(BasisFuncND("bf_"+std::to_string(id),_n_dim,_dims));
			// Updates to the basis functions
			_dbfs.push_back(GridND("dbf_"+std::to_string(id),_n_dim,_dims));
		};

		// Initialize solution trajs
		if (DIAG_INIT) {
			std::cout << "soln traj..." << std::flush;
		};
		_soln_traj = new double*[n_t];
		for (int t=0; t<n_t; t++) {
			_soln_traj[t] = new double[_n_dim];
		};

		// Variational problem solution trajectories
		if (DIAG_INIT) {
			std::cout << "var problem..." << std::flush;
		};
		for (int t=0; t<n_t; t++) 
		{
			_var_traj.push_back(std::vector<std::vector<GridND>>());
			for (int d_num=0; d_num<_n_dim; d_num++) 
			{
				_var_traj[t].push_back(std::vector<GridND>());
				for (int d_denom=0; d_denom<_n_dim; d_denom++) 
				{
					_var_traj[t][d_num].push_back(GridND("var_t_"+std::to_string(t)+"_vary_"+std::to_string(d_num)+"_wrt_"+std::to_string(d_denom),_n_dim,_dims));
				};
			};
		};

		// Species
		if (DIAG_INIT) {
			std::cout << "species..." << std::flush;
		};
		for (auto sp: species) {
			// Add the species
			_species.push_back(Species(sp));
			// Link ptrs
			_species.back()._soln_traj_ptr = &_soln_traj;
			_species.back()._t_opt_ptr = &_t_opt;

			// Add to the lattice
			_latt.add_species(&(_species.back()));
		};

		// Go through the dimensions and link them to the species
		if (DIAG_INIT) {
			std::cout << "linking dim->species..." << std::flush;
		};
		for (int d=0; d<_n_dim; d++)
		{
			if (_dims[d].type == H) {
				_dims[d].sp1 = _find_species_by_name(_dims[d].species1);
			} else if (_dims[d].type == J) {
				_dims[d].sp1 = _find_species_by_name(_dims[d].species1);
				_dims[d].sp2 = _find_species_by_name(_dims[d].species2);
			};
		};

		// Go through the species and link them to the params
		if (DIAG_INIT) {
			std::cout << "linking species->params..." << std::flush;
		};
		for (auto its1=_species.begin(); its1!=_species.end(); its1++)
		{
			its1->_h_index = _find_param_by_species(&*its1);
			for (auto its2=_species.begin(); its2!=_species.end(); its2++)
			{
				its1->_j_index[&*its2] = _find_param_by_species(&*its1,&*its2);
			};
		};

		// Optimization params
		this->_n_batch = batch_size;
		this->_n_annealing = n_annealing;
		this->_box_length = box_length;
		this->_dopt = dopt;
		this->_n_opt = n_opt; 

		// Allocate data structures to solve variational term
		if (DIAG_INIT) {
			std::cout << "allocating structures for var problem..." << std::flush;
		};
		_bfs_derivs_var = new double*[_n_dim];
		for (int d=0; d<_n_dim; d++) {
			_bfs_derivs_var[d] = new double[_n_dim];
		};
		_var_idxs = new int[_n_dim];

		if (DIAG_INIT) {
			std::cout << "done." << std::endl;
		};
	};
	OptProblem::~OptProblem() {
		for (int t=0; t<_time.n; t++) {
			safeDelArr(_soln_traj[t]);
		};
		delete[] _soln_traj;

		safeDelArr(_var_idxs);

		for (int d=0; d<_n_dim; d++) {
			safeDelArr(_bfs_derivs_var[d]);
		};
		delete[] _bfs_derivs_var;
	};

	/********************
	Set properties
	********************/

	void OptProblem::add_fname(std::string f) {
		_fnames.push_back(f);
	};

	/********************
	Solve for solution traj
	********************/

	void OptProblem::solve_traj()
	{
		// Initialize
		for (int d=0; d<_n_dim; d++)
		{
			_soln_traj[0][d] = _init[d];
		};

		// The max time the integration has gotten to
		_n_t_soln = 0;

		// The new point to add
		double pnew;

		// Go through all times
		for (int t=0; t<_time.n-1; t++) {
			// All dims
			for (int d=0; d<_n_dim; d++) {
				pnew = _soln_traj[t][d] + _time.delta * _bfs[d].get_val(_soln_traj[t]);
				// Check that the new point has not gone out of bounds
				if (pnew > _dims[d].max || pnew < _dims[d].min) {
					std::cout << "WARNING: traj has gone out of bounds. dim " << _dims[d].name << " point " << pnew << " out of [" << _dims[d].min << "," << _dims[d].max << "]." << std::endl;
					// Trajectory has gone out of the bounds of the grid
					return;
				};
				// Add the point
				_soln_traj[t+1][d] = pnew;
			};
			// All dims are done for this time
			_n_t_soln++;
		};
	};

	/********************
	Solve for var traj
	********************/

	void OptProblem::solve_var_traj()
	{
		// Initial = all zero

		if (DIAG_VAR) { std::cout << "Init var..." << std::flush; };

		for (int d_num=0; d_num<_n_dim; d_num++) {
			for (int d_denom=0; d_denom<_n_dim; d_denom++) {
				_var_traj[0][d_num][d_denom].zero();
			};
		};

		if (DIAG_VAR) { std::cout << "ok." << std::endl; };

		// To specify which dimension to differentiate
		bool derivs[_n_dim];
		std::fill_n(derivs,_n_dim,false); // all false

		// Go through all times, as long as the soln traj has not gone out of bounds
		for (_var_time=0; _var_time<_n_t_soln-1; _var_time++) {

			// Fill out the derivatives of the basis function at the soln point corresponding to this time

			if (DIAG_VAR) { std::cout << "Fill out derivs at time " << _var_time << "..." << std::flush; };

			for (int d_deriv=0; d_deriv<_n_dim; d_deriv++) { // Deriv dimension
				derivs[d_deriv] = true; // diff this dimension
				for (int d_bf=0; d_bf<_n_dim; d_bf++) { // BF dimension
					_bfs_derivs_var[d_bf][d_deriv] = _bfs[d_bf].get_deriv(_soln_traj[_var_time],derivs);
				};
				derivs[d_deriv] = false; // reset
			};

			if (DIAG_VAR) { std::cout << "ok." << std::endl; };

			// Iterate over all variational terms

			if (DIAG_VAR) { std::cout << "Iterating over grid..." << std::endl; };

			for (_var_nu=0; _var_nu<_n_dim; _var_nu++) {
				for (_var_f=0; _var_f<_n_dim; _var_f++) {					
					// Recursion to iterate over all grid points
					_iterate_var(0);
				};
			};

			if (DIAG_VAR) { std::cout << "ok." << std::endl; };

		};
	};

	// Recursive function to iterate over all grid points for the variational term
	void OptProblem::_iterate_var(int dim)
	{
		// Iterate over all idxs in dimension dim
		for(_var_idxs[dim] = 0; _var_idxs[dim] < _dims[dim].n; ++_var_idxs[dim]) 
		{
			if (dim != _n_dim-1)
			{
				// Inception - need another for loop
				_iterate_var(dim+1);
			} else {
				// Do something!

				// Val at prev time
				if (DIAG_VAR) { 
					std::cout << "   Doing point: ";
					for (int d=0; d<_n_dim; d++) {
						std::cout << _dims[d].grid[_var_idxs[d]] << " ";
					};
					std::cout << "..." << std::flush; 
				};

				double new_val = _var_traj[_var_time][_var_nu][_var_f].get(_var_idxs);

				// Delta source
				if (_var_nu == _var_f) {
					new_val += _time.delta * _var_delta();
				};

				// Derivatives
				// We are calculating:
				// (delta _var_nu) / (delta _var_f)
				// These terms are: 
				// sum_l (delta nu_l) / (delta _var_f) * (deriv of F corresponding to _var_nu ) / (d nu_l)
				// So: _n_dim terms in the sum!
				for (int nu_l=0; nu_l<_n_dim; nu_l++) {
					new_val += _time.delta * _var_traj[_var_time][nu_l][_var_f].get(_var_idxs) * _bfs_derivs_var[_var_nu][nu_l];
				};

				if (DIAG_VAR) { std::cout << " val " << new_val << std::flush; };

				// Set val at new time
				_var_traj[_var_time+1][_var_nu][_var_f].set(_var_idxs,new_val);

				if (DIAG_VAR) { std::cout << " ok." << std::endl; };
			};
		};
	};	


	/********************
	Write
	********************/

	void OptProblem::write_soln_traj(std::string fname) {
		std::ofstream f;
		f.open (fname);
		for (int i=0; i<_n_t_soln; i++) {
			f << _time.grid[i] << " ";
			for (int d=0; d<_n_dim; d++) {
				f << _soln_traj[i][d];
				if (d != _n_dim-1) { f << " "; };
			};
			if (i != _n_t_soln-1) { f << "\n";};
		};
		f.close();
	};
	void OptProblem::write_var_traj(std::string fname) {
		// Length of the grid
		int val_len = 1;
		for (auto v: _dims) { val_len *= v.n; };

		std::ofstream f;
		f.open (fname);
		for (int t=0; t<_n_t_soln; t++) {
			// Loop over grid
			for (int i=0; i<val_len; i++) {
				f << _time.grid[t] << " ";
				// Loop over num/denom and write value
				for (int d_num=0; d_num<_n_dim; d_num++) {
					for (int d_denom=0; d_denom<_n_dim; d_denom++) {
						f << _var_traj[t][d_num][d_denom].get(i) << " ";
					};
				};
				f << "\n";
			};
		};
		f.close();
	};
	void OptProblem::write_bfs(std::string fname) {
		// Length of the grid
		int val_len = 1;
		for (auto v: _dims) { val_len *= v.n; };

		std::ofstream f;
		f.open (fname);
		// Loop over grid
		for (int i=0; i<val_len; i++) {
			// Loop over basis funcs
			for (int d=0; d<_n_dim; d++) {
				f << _bfs[d].get(i) << " ";
			};
			f << "\n";
		};
		f.close();
	};
	void OptProblem::write_bfs(std::string dir,int idx) {
		for (int d=0; d<_n_dim; d++) 
		{
			_bfs[d].write_to_file(dir,idx);
		};
	};

	void OptProblem::write_grid(std::string fname) {
		// Length of the grid
		int val_len = 1;
		for (auto v: _dims) { val_len *= v.n; };

		// Idxs
		int *idxs = new int[_n_dim];

		std::ofstream f;
		f.open (fname);
		// Loop over grid
		for (int i=0; i<val_len; i++) {
			// Get idxs
			idxs = _var_traj[0][0][0].get_idxs(i);
			// Write coords
			for (int d=0; d<_n_dim; d++) {
				f << _dims[d].grid[idxs[d]] << " ";
			};
			if (i != val_len-1) { f << "\n"; };
		};
		f.close();
		
		delete[] idxs;	
	};

	void OptProblem::write_moments(std::string fname, bool append) {
		std::ofstream f;
		if (append) {
			f.open(fname,std::ofstream::out | std::ofstream::app);
		} else {
			f.open(fname);
		};
		// Write time
		f << _t_opt << " ";
		// Loop over dims
		for (int d=0; d<_n_dim; d++) {
			f << _dims[d].awake / _n_batch << " ";
		};
		for (int d=0; d<_n_dim; d++) {
			f << _dims[d].asleep / _n_batch << " ";
		};
		f << "\n";
		f.close();
	};

	/********************
	Run main optimization loop
	********************/

	void OptProblem::solve(bool verbose) {

		// Write the grid of points
		if (DIAG_SOLVE) { std::cout << "Writing grid..." << std::flush; };
		write_grid("data/grid.txt");
		if (DIAG_SOLVE) { std::cout << "ok." << std::endl; };

		// Filename to read
		std::string fname;

		// Multiplicative factor
		double m_factor;

		// Iterate over optimization steps
		for (int i_opt=0; i_opt<_n_opt; i_opt++)
		{
			std::cout << "Opt step " << i_opt << " / " << _n_opt-1 << std::endl;

			// Write the current basis funcs

			if (DIAG_SOLVE) { std::cout << "Writing F..." << std::flush; };

			write_bfs("data/F/"+pad_str(i_opt,4)+".txt");
			
			if (DIAG_SOLVE) { std::cout << "ok." << std::endl; };

			/*****
			Step 1 - Solve the current trajectory
			*****/

			if (DIAG_SOLVE) { std::cout << "Solving..." << std::flush; };

			solve_traj();
			
			if (DIAG_SOLVE) { std::cout << "ok." << std::endl; };

			// Write the trajectory

			if (DIAG_SOLVE) { std::cout << "Writing soln..." << std::flush; };	
			write_soln_traj("data/soln_traj/"+pad_str(i_opt,4)+".txt");

			if (DIAG_SOLVE) { std::cout << "ok." << std::endl; };

			/*****
			Step 2 - Solve the variational problem traj
			*****/

			if (DIAG_SOLVE) { std::cout << "Solving var..." << std::flush; };

			solve_var_traj();
			
			if (DIAG_SOLVE) { std::cout << "ok." << std::endl; };

			// Write the nu trajectory

			if (DIAG_SOLVE) { std::cout << "Writing var..." << std::flush; };
			
			write_var_traj("data/var_traj/"+pad_str(i_opt,4)+".txt");
			
			if (DIAG_SOLVE) { std::cout << "ok." << std::endl; };

			/*****
			Step 3
			*****/

			// Pick a random batch

			if (DIAG_SOLVE) { std::cout << "Got batch..." << std::flush; };

			std::vector<std::string> fnames, fnames_possible=_fnames;
			std::vector<std::string>::iterator itf;
			for (int i_batch=0; i_batch<_n_batch; i_batch++) {
				itf = fnames_possible.begin();
				std::advance(itf,randI(0,fnames_possible.size()-1));
				fnames.push_back(*itf);
				fnames_possible.erase(itf);
			};

			if (DIAG_SOLVE) { std::cout << "ok." << std::endl; };

			/*****
			Step 4
			*****/

			// Reset update
			for (int d=0; d<_n_dim; d++) {
				_dbfs[d].zero();
			};

			/*****
			Step 5
			*****/

			// Go through all times
			// Use the class variable _t_opt to iterate
			for (_t_opt=0; _t_opt < _n_t_soln; _t_opt++)
			{
				if (verbose) {
					std::cout << "time: " << _t_opt << std::flush;
				};

				/*****
				Step 5.1
				*****/

				// Reset the moments to zero
				for (int d=0; d<_n_dim; d++) {
					_dims[d].awake = 0.0;
					_dims[d].asleep = 0.0;
				};

				/*****
				Step 5.2
				*****/

				// Go through all samples in the batch
				for (int i_batch=0; i_batch<_n_batch; i_batch++) 
				{
					if (verbose) {
						std::cout << "." << std::flush;
					};

					/*****
					Step 5.2.1
					*****/

					// Read in batch at this timestep
					fname = fnames[i_batch] + pad_str(_t_opt,4) + ".txt";
					_latt.read_from_file(fname);

					/*****
					Step 5.2.2
					*****/

					// Record the awake moments
					for (int d=0; d<_n_dim; d++) {
						_dims[d].append_moments_from_latt(Dim::AWAKE);
					};

					/*****
					Step 5.2.3
					*****/

					// Anneal
					_latt.anneal(_n_annealing);

					/*****
					Step 5.2.4
					*****/

					// Record the asleep moments
					for (int d=0; d<_n_dim; d++) {
						_dims[d].append_moments_from_latt(Dim::ASLEEP);
					};
				};

				// Write the moments
				if (_t_opt == 0) {
					write_moments("data/moments/"+pad_str(i_opt,4)+".txt",false);
				} else {
					write_moments("data/moments/"+pad_str(i_opt,4)+".txt",true);
				};

				/*****
				Step 5.3
				*****/

				// Write the moments into the _var_traj
				// This invalidates the _var_traj but is more efficient!
				// Update at this time looks like:
				// sum_l ( delta S / delta nu_l ) * ( delta nu_l / delta F )
				// Multiply by the first term
				for (int l=0; l<_n_dim; l++) { // loop numerator delta nu
					m_factor = -1.0 * _dopt * _time.delta * (_dims[l].asleep - _dims[l].awake) / _n_batch;
					for (int d_bf=0; d_bf<_n_dim; d_bf++) { // loop denom delta F
						_var_traj[_t_opt][l][d_bf].multiply_all(m_factor);
					};
				};

				/*****
				Step 5.4
				*****/

				// Gather the update

				// Loop over basis funcs to update
				for (int d_bf=0; d_bf<_n_dim; d_bf++)
				{
					// Update at this time looks like:
					// sum_l ( delta S / delta nu_l ) * ( delta nu_l / delta F )
					// These are now in _var_traj
					// Do sum over l
					for (int l=0; l<_n_dim; l++) {
						_dbfs[d_bf].increment(_var_traj[_t_opt][l][d_bf]);
					};
				};

				if (verbose) {
					std::cout << "done." << std::endl;
				};
			};

			/*****
			Step 6
			*****/

			// Update the basis funcs
			for (int d_bf=0; d_bf<_n_dim; d_bf++)
			{
				_bfs[d_bf].update(_dbfs[d_bf]);
			};
		};

		// Write the final basis funcs once more
		write_bfs("data/F/"+pad_str(_n_opt,4)+".txt");
	};

	/****************************************
	OptProblem - PRIVATE
	****************************************/

	/********************
	Delta function at some time index, nu value
	********************/

	double OptProblem::_var_delta()
	{
		double r=1;
		for (int d=0; d<_n_dim; d++) {
			r *= exp(-pow(_soln_traj[_var_time][d]-_dims[d].grid[_var_idxs[d]],2)/(2.*0.1))/sqrt(2.*M_PI*0.1);
		};
		return r;
	};

	/********************
	Find a...
	********************/

	Species* OptProblem::_find_species_by_name(std::string name) {
		// Go through the species
		auto it = this->_species.begin();
		while (it != this->_species.end()) {
			if (it->name == name) {
				return &(*it);
			};
			it++;
		};
		return nullptr;
	};

	int OptProblem::_find_param_by_species(Species *sp) {
		for (int d=0; d<_n_dim; d++) {
			if (_dims[d].type == H) {
				if (_dims[d].sp1 == sp) {
					return d;
				};
			};
		};
		return -1;
	};

	int OptProblem::_find_param_by_species(Species *sp1, Species *sp2) {
		for (int d=0; d<_n_dim; d++) {
			if (_dims[d].type == J) {
				if ( (_dims[d].sp1 == sp1 && _dims[d].sp2 == sp2) || (_dims[d].sp1 == sp2 && _dims[d].sp2 == sp1) ) {
					return d;
				};
			};
		};
		return -1;
	};
};


