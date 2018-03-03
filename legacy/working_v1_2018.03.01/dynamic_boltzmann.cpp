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

		// Update the dimensions for each basis function
		if (DIAG_INIT) {
			std::cout << "updating basis funcs of dims..." << std::flush;
		};
		std::vector<int> idxs_all(_n_dim);
		std::iota(idxs_all.begin(), idxs_all.end(), 0);
		for (int d=0; d<_n_dim; d++) {
			// By default, use all dims for the basis func corresponding to this dim
			if (_dims[d].bf_use_all_dim) {
				_dims[d].bf_n_dim = _n_dim;
				_dims[d].bf_dim_idxs = idxs_all;
			} else { // Otherwise, use the given dims
				_dims[d].bf_dim_idxs.clear();
				for (auto dim_name: _dims[d].bf_dim_names)
				{
					// Search for the idx of the dimension with this name
					int dsearch=0; bool found = false;
					while (dsearch < _n_dim && !found) {
						if (_dims[dsearch].name==dim_name) {
							_dims[d].bf_dim_idxs.push_back(dsearch);
							found = true;
						} else {
							dsearch++;
						};
					};
					if (!found) {
						std::cerr << "ERROR! Dimension with name: " << dim_name << " was not found." << std::endl;
						exit(EXIT_FAILURE);
					};
				};
			};
		};

		// Initial points
		_init = init;

		// Initialize basis funcs and updates
		if (DIAG_INIT) {
			std::cout << "basis funcs..." << std::flush;
		};
		std::vector<Dim> bf_dims;
		for (int d=0; d<_n_dim; d++) {
			// Set up the dimensions
			bf_dims.clear();
			for (auto dim_bf: _dims[d].bf_dim_idxs) {
				bf_dims.push_back(_dims[dim_bf]);
			};
			// Make basis func
			_bfs.push_back(BasisFuncND("F_"+_dims[d].name,_dims[d].bf_n_dim,bf_dims));
			_dbfs.push_back(GridND("delta_F_"+_dims[d].name,_dims[d].bf_n_dim,bf_dims));
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
					// Set up the dimensions
					bf_dims.clear();
					for (auto dim_bf: _dims[d_denom].bf_dim_idxs) {
						bf_dims.push_back(_dims[dim_bf]);
					};
					_var_traj[t][d_num].push_back(GridND("var_t_"+std::to_string(t)+"_vary_"+_dims[d_num].name+"_wrt_"+_dims[d_denom].name,_dims[d_denom].bf_n_dim,bf_dims));
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
		_var_bf_derivs = new double*[_n_dim];
		for (int d=0; d<_n_dim; d++) {
			_var_bf_derivs[d] = new double[_dims[d].bf_n_dim];
		};
		// Do not init _var_idxs - its size changes

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
			safeDelArr(_var_bf_derivs[d]);
		};
		delete[] _var_bf_derivs;
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
		_n_t_soln = 1;

		// The new point to add
		double pnew;

		// Current point
		double **pcurr;
		pcurr = new double*[_n_dim];
		for (int d=0; d<_n_dim; d++) {
			pcurr[d] = new double[_dims[d].bf_n_dim];
		};

		// Go through all times
		for (int t=0; t<_time.n-1; t++) {

			if (DIAG_SOLVE_TRAJ) { std::cout << t << " " << std::flush; };

			// All dims
			for (int d=0; d<_n_dim; d++) {

				// Point for this basis func
				for (int i=0; i<_dims[d].bf_n_dim; i++)  {
					pcurr[d][i] = _soln_traj[t][_dims[d].bf_dim_idxs[i]];
				};

				// New point
				pnew = _soln_traj[t][d] + _time.delta * _bfs[d].get_val(pcurr[d]);
				
				// Check that the new point has not gone out of bounds
				if (pnew > _dims[d].max || pnew < _dims[d].min) {
					std::cout << "WARNING: traj has gone out of bounds. dim " << _dims[d].name << " point " << pnew << " out of [" << _dims[d].min << "," << _dims[d].max << "]." << std::endl;
					// Trajectory has gone out of the bounds of the grid
					// Free
					for (int d2=0; d2<_n_dim; d2++) {
						safeDelArr(pcurr[d2]);
					};
					delete[] pcurr;

					return;
				};
				// Add the point
				_soln_traj[t+1][d] = pnew;

				if (DIAG_SOLVE_TRAJ) { std::cout << pnew << " " << std::flush; };
			};

			if (DIAG_SOLVE_TRAJ) { std::cout << std::endl; };

			// All dims are done for this time
			_n_t_soln++;
		};

		// Free
		for (int d2=0; d2<_n_dim; d2++) {
			safeDelArr(pcurr[d2]);
		};
		delete[] pcurr;
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
		bool *derivs;

		// Current point
		double **pcurr;
		pcurr = new double*[_n_dim];
		for (int d=0; d<_n_dim; d++) {
			pcurr[d] = new double[_dims[d].bf_n_dim];
		};

		// Go through all times, as long as the soln traj has not gone out of bounds
		for (_var_time=0; _var_time<_n_t_soln-1; _var_time++) {

			// Fill out the derivatives of the basis function at the soln point corresponding to this time

			if (DIAG_VAR) { std::cout << "Fill out derivs at time " << _var_time << "..." << std::flush; };

			// Go through basis funcs to differentiate
			for (int d_bf=0; d_bf<_n_dim; d_bf++) {

				// Point for this basis func
				for (int i=0; i<_dims[d_bf].bf_n_dim; i++)  {
					pcurr[d_bf][i] = _soln_traj[_var_time][_dims[d_bf].bf_dim_idxs[i]];
				};

				// The deriv to take
				derivs = new bool[_dims[d_bf].bf_n_dim];
				std::fill_n(derivs,_dims[d_bf].bf_n_dim,false); // all false

				// Go through nu to differentiate with respect to
				for (int d_nu=0; d_nu<_dims[d_bf].bf_n_dim; d_nu++) {
					derivs[d_nu] = true; // diff this dimension
					_var_bf_derivs[d_bf][d_nu] = _bfs[d_bf].get_deriv(pcurr[d_bf],derivs);
					derivs[d_nu] = false; // reset
				};
				// Free
				safeDelArr(derivs);
			};

			if (DIAG_VAR) { std::cout << "ok." << std::endl; };

			// Iterate over all variational terms

			if (DIAG_VAR) { std::cout << "Iterating over grid..." << std::endl; };

			for (_var_nu=0; _var_nu<_n_dim; _var_nu++) {
				for (_var_f=0; _var_f<_n_dim; _var_f++) {
					// Init the idxs for the grid of the basis function we are varying with respect to
					_var_f_n_dim = _dims[_var_f].bf_n_dim;
					_var_nu_n_dim = _dims[_var_nu].bf_n_dim;
					_var_idxs = new int[_var_f_n_dim];
					_var_f_dim_idxs = _dims[_var_f].bf_dim_idxs;
					_var_nu_dim_idxs = _dims[_var_nu].bf_dim_idxs;

					// Recursion to iterate over all grid points
					_iterate_var(0);

					// Free the idxs
					safeDelArr(_var_idxs);
				};
			};

			if (DIAG_VAR) { std::cout << "ok." << std::endl; };

		};

		// Free
		for (int d2=0; d2<_n_dim; d2++) {
			safeDelArr(pcurr[d2]);
		};
		delete[] pcurr;
	};

	// Recursive function to iterate over all grid points for the variational term
	void OptProblem::_iterate_var(int dim)
	{
		// Iterate over all idxs in dimension dim
		for(_var_idxs[dim] = 0; _var_idxs[dim] < _dims[_var_f_dim_idxs[dim]].n; ++_var_idxs[dim]) 
		{
			if (dim != _var_f_n_dim-1)
			{
				// Inception - need another for loop
				_iterate_var(dim+1);
			} else {
				// Do something!

				if (DIAG_VAR) { 
					std::cout << "   Doing point: ";
					for (int d=0; d<_var_f_n_dim; d++) {
						std::cout << _dims[_var_f_dim_idxs[d]].grid[_var_idxs[d]] << " ";
					};
					std::cout << "..." << std::flush; 
				};

				// Val at prev time
				double new_val = _var_traj[_var_time][_var_nu][_var_f].get(_var_idxs);
				//if (new_val > 1e10) { std::cerr << "prev" << std::endl; exit(EXIT_FAILURE); };

				// Delta source
				if (_var_nu == _var_f) {
					new_val += _time.delta * _var_delta();
				};
				//if (new_val > 1e10) { std::cerr << "delta" << std::endl; exit(EXIT_FAILURE); };

				// Derivatives
				// We are calculating:
				// (delta _var_nu) / (delta _var_f)
				// These terms are: 
				// sum_l (delta nu_l) / (delta _var_f) * (deriv of F corresponding to _var_nu ) / (d nu_l)
				// l runs over the number of basis funcs for this dim
				// So: _var_f_n_dim terms in the sum!
				for (int nu_l_idx=0; nu_l_idx<_var_nu_n_dim; nu_l_idx++) {
					new_val += _time.delta * _var_traj[_var_time][_var_nu_dim_idxs[nu_l_idx]][_var_f].get(_var_idxs) * _var_bf_derivs[_var_nu][nu_l_idx];
					// if (new_val > 1e10) { std::cerr << "deriv " << _var_time << " " << _var_nu << " " << _var_f << " " << nu_l_idx << " " << _var_nu_dim_idxs[nu_l_idx] << " " << _var_traj[_var_time][_var_nu_dim_idxs[nu_l_idx]][_var_f].get(_var_idxs) << " " << _var_bf_derivs[_var_nu][nu_l_idx] << std::endl; exit(EXIT_FAILURE); };
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
			for (int d=0; d<_n_dim; d++) {
				f << _soln_traj[i][d];
				if (d != _n_dim-1) { f << " "; };
			};
			if (i != _n_t_soln-1) { f << "\n";};
		};
		f.close();
	};
	void OptProblem::write_var_traj(std::string fname, int idx_num, int idx_denom) {
		// Length of the grid
		int val_len = 1;
		for (auto i: _dims[idx_denom].bf_dim_idxs) { val_len *= _dims[i].n; };

		std::ofstream f;
		f.open (fname);
		// Loop over time
		for (int t=0; t<_n_t_soln; t++) {
			// Loop over grid
			for (int i=0; i<val_len; i++) {
				f << _var_traj[t][idx_num][idx_denom].get(i) << "\n";
			};
		};
		f.close();
	};
	void OptProblem::write_bf(std::string fname, int idx) {
		// Length of the grid
		int val_len = 1;
		for (auto i: _dims[idx].bf_dim_idxs) { val_len *= _dims[i].n; };

		std::ofstream f;
		f.open (fname);
		// Loop over grid
		for (int i=0; i<val_len; i++) {
			f << _bfs[idx].get(i);
			if (i != val_len-1) { f << "\n"; };
		};
		f.close();
	};
	void OptProblem::write_bf_grid(std::string fname, int idx) {
		// Length of the grid
		int val_len = 1;
		for (auto i: _dims[idx].bf_dim_idxs) { val_len *= _dims[i].n; };

		// Idxs
		int *idxs = new int[_dims[idx].bf_n_dim];

		std::ofstream f;
		f.open (fname);
		// Loop over grid
		for (int i=0; i<val_len; i++) {
			// Get idxs
			idxs = _bfs[idx].get_idxs(i);
			// Write coords
			for (int d=0; d<_dims[idx].bf_n_dim; d++) {
				f << _dims[_dims[idx].bf_dim_idxs[d]].grid[idxs[d]] << " ";
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

	void OptProblem::write_t_grid(std::string fname) {
		std::ofstream f;
		f.open(fname);
		for (int t=0; t<_time.n; t++) {
			f << _time.grid[t];
			if (t != _time.n-1) { f << "\n"; };
		};
		f.close();
	};

	/********************
	Run main optimization loop
	********************/

	void OptProblem::solve(bool verbose) {

		// Write the grid of points
		if (DIAG_SOLVE) { std::cout << "Writing grid..." << std::flush; };
		for (int i=0; i<_n_dim; i++) {
			write_bf_grid("data/grid_"+_bfs[i].name+".txt",i);
		};
		write_t_grid("data/grid_time.txt");
		if (DIAG_SOLVE) { std::cout << "ok." << std::endl; };

		// Multiplicative factor
		double m_factor;

		// Iterate over optimization steps
		for (int i_opt=0; i_opt<_n_opt; i_opt++)
		{
			std::cout << "Opt step " << i_opt << " / " << _n_opt-1 << std::endl;

			// Write the current basis funcs

			if (DIAG_SOLVE) { std::cout << "Writing F..." << std::flush; };

			for (int i=0; i<_n_dim; i++) {
				write_bf("data/F/"+_bfs[i].name+"_"+pad_str(i_opt,4)+".txt",i);
			};
			
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

			// Write the trajectory

			if (DIAG_SOLVE) { std::cout << "Writing var..." << std::flush; };
			
			/*
			for (int i_num=0; i_num<_n_dim; i_num++) {
				for (int i_denom=0; i_denom<_n_dim; i_denom++) {
					write_var_traj("data/var_traj/"+_dims[i_num].name+"_"+_dims[i_denom].name+"_"+pad_str(i_opt,4)+".txt",i_num,i_denom);
				};
			};
			*/
			
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
					_latt.read_from_file(fnames[i_batch] + pad_str(_t_opt,4) + ".txt");

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
		for (int i=0; i<_n_dim; i++) {
			write_bf("data/F/"+_bfs[i].name+"_"+pad_str(_n_opt,4)+".txt",i);
		};
	};

	/********************
	Solve for varying initial conditions
	********************/

	void OptProblem::solve_varying_ic(bool verbose) {

		// Write the grid of points
		if (DIAG_SOLVE) { std::cout << "Writing grid..." << std::flush; };
		for (int i=0; i<_n_dim; i++) {
			write_bf_grid("data/grid_"+_bfs[i].name+".txt",i);
		};
		write_t_grid("data/grid_time.txt");
		if (DIAG_SOLVE) { std::cout << "ok." << std::endl; };

		// Multiplicative factor
		double m_factor;

		// To write out the moments
		double **moms_awake_st;
		double **moms_asleep_st;
		moms_awake_st = new double*[_time.n];
		moms_asleep_st = new double*[_time.n];
		for (int t=0; t<_time.n; t++) {
			moms_awake_st[t] = new double[_n_dim];
			moms_asleep_st[t] = new double[_n_dim];
		};

		// Iterate over optimization steps
		for (int i_opt=0; i_opt<_n_opt; i_opt++)
		{
			std::cout << "Opt step " << i_opt << " / " << _n_opt-1 << std::endl;

			// Write the current basis funcs

			if (DIAG_SOLVE) { std::cout << "Writing F..." << std::flush; };

			for (int i=0; i<_n_dim; i++) {
				write_bf("data/F/"+_bfs[i].name+"_"+pad_str(i_opt,4)+".txt",i);
			};
			
			if (DIAG_SOLVE) { std::cout << "ok." << std::endl; };

			/*****
			Step 1
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
			Step 2
			*****/

			// Reset update
			for (int d=0; d<_n_dim; d++) {
				_dbfs[d].zero();
			};

			/*****
			Step 3
			*****/

			// Go through all samples in the batch
			for (int i_batch=0; i_batch<_n_batch; i_batch++) 
			{
				if (verbose) {
					std::cout << "batch: " << i_batch << " / " << _n_batch-1 << std::flush;
				};

				/*****
				Step 3.1 - Read in the initial conditions
				*****/

				if (DIAG_SOLVE) { std::cout << "Reading IC..." << std::flush; };

				read_ic(fnames[i_batch] + "init.txt");

				if (DIAG_SOLVE) { std::cout << "ok." << std::endl; };

				/*****
				Step 3.2 - Solve the current trajectory
				*****/

				if (DIAG_SOLVE) { std::cout << "Solving..." << std::flush; };

				solve_traj();
				
				if (DIAG_SOLVE) { std::cout << "ok." << std::endl; };

				// Write the trajectory

				if (DIAG_SOLVE) { std::cout << "Writing soln..." << std::flush; };	

				write_soln_traj("data/soln_traj/"+pad_str(i_opt,4)+"_"+pad_str(i_batch,3)+".txt");

				if (DIAG_SOLVE) { std::cout << "ok." << std::endl; };

				/*****
				Step 3.3 - Solve the variational problem traj
				*****/

				if (DIAG_SOLVE) { std::cout << "Solving var..." << std::flush; };

				solve_var_traj();
				
				if (DIAG_SOLVE) { std::cout << "ok." << std::endl; };

				// Write the trajectory

				if (DIAG_SOLVE) { std::cout << "Writing var..." << std::flush; };
				
				/*
				for (int i_num=0; i_num<_n_dim; i_num++) {
					for (int i_denom=0; i_denom<_n_dim; i_denom++) {
						write_var_traj("data/var_traj/"+_dims[i_num].name+"_"+_dims[i_denom].name+"_"+pad_str(i_opt,4)+".txt",i_num,i_denom);
					};
				};
				*/

				if (DIAG_SOLVE) { std::cout << "ok." << std::endl; };

				/*****
				Step 3.4
				*****/

				// Go through all times
				// Use the class variable _t_opt to iterate
				for (_t_opt=0; _t_opt < _n_t_soln; _t_opt++)
				{
					if (verbose) {
						std::cout << "." << std::flush;
					};

					/*****
					Step 3.4.1
					*****/

					// Reset the moments to zero
					for (int d=0; d<_n_dim; d++) {
						_dims[d].awake = 0.0;
						_dims[d].asleep = 0.0;
					};

					/*****
					Step 3.4.2
					*****/

					// Read in batch at this timestep
					_latt.read_from_file(fnames[i_batch] + pad_str(_t_opt,4) + ".txt");

					/*****
					Step 3.4.3
					*****/

					// Record the awake moments
					for (int d=0; d<_n_dim; d++) {
						_dims[d].append_moments_from_latt(Dim::AWAKE);
						moms_awake_st[_t_opt][d] = _dims[d].awake;
					};

					/*****
					Step 3.4.4
					*****/

					// Anneal
					_latt.anneal(_n_annealing);

					/*****
					Step 3.4.5
					*****/

					// Record the asleep moments
					for (int d=0; d<_n_dim; d++) {
						_dims[d].append_moments_from_latt(Dim::ASLEEP);
						moms_asleep_st[_t_opt][d] = _dims[d].asleep;
					};

					/*****
					Step 3.4.6
					*****/

					// Write the moments into the _var_traj
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
					Step 3.4.7
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

					/*****
					3.4.8
					*****/

					// Reset the var traj
					for (int l=0; l<_n_dim; l++) { // loop numerator delta nu
						m_factor = -1.0 * _dopt * _time.delta * (_dims[l].asleep - _dims[l].awake) / _n_batch;
						for (int d_bf=0; d_bf<_n_dim; d_bf++) { // loop denom delta F
							_var_traj[_t_opt][l][d_bf].multiply_all(1.0/m_factor);
						};
					};
				};

				// Write the moments for this batch
				std::ofstream f;
				f.open("data/moments/"+pad_str(i_opt,4)+"_"+pad_str(i_batch,3)+".txt");
				// Loop over time
				for (int t=0; t<_n_t_soln; t++) {
					// Loop over dims
					for (int d=0; d<_n_dim; d++) {
						f << moms_awake_st[t][d] << " ";
					};
					for (int d=0; d<_n_dim; d++) {
						f << moms_asleep_st[t][d];
						if (d != _n_dim-1) { f << " "; };
					};
					if (t != _n_t_soln-1) { f << "\n"; };
				};
				f.close();

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
		for (int i=0; i<_n_dim; i++) {
			write_bf("data/F/"+_bfs[i].name+"_"+pad_str(_n_opt,4)+".txt",i);
		};

		// Clean up
		for (int t=0; t<_time.n; t++) {
			delete[] moms_awake_st[t];
			delete[] moms_asleep_st[t];
		};
		delete[] moms_awake_st;
		delete[] moms_asleep_st;
	};

	// Read in initial conditions
	void OptProblem::read_ic(std::string fname) {

		std::ifstream f;
		f.open(fname);
		char frag[100]; // fragments of the line
		int i_frag=0;
		int i_line=0;
		std::string s1="",s2="";
		std::string val="";

		std::map<std::string,double> h0; 
		std::map<std::string,std::map<std::string,double>> j0;

		if (f.is_open()) { // make sure we found it
			while (!f.eof()) {
			    f >> frag;
			    if (i_line < 3) {
			    	if (i_frag==0) { 
			    		s1 += frag; i_frag++;
			    	} else if (i_frag==1) {
			    		val += frag;
			    		// Add to map
			    		h0[s1] = std::stod(val);
			    		// Reset
			    		i_line++;
			    		i_frag=0;
			    		s1="";s2="";val="";
			    	};
			    } else {
			    	if (i_frag==0) { 
			    		s1 += frag; i_frag++;
			    	} else if (i_frag==1) { 
			    		s2 += frag; i_frag++;
			    	} else if (i_frag==2) {
			    		val += frag;
			    		// Add to map
			    		j0[s1][s2] = std::stod(val);
			    		// Reset
			    		i_line++;
			    		i_frag=0;
			    		s1="";s2="";val="";
			    	};
			    };
			};
		};

		// Translate the maps into the _init vector
		_init.clear();
		bool swtch = false;
		for (int d=0; d<_n_dim; d++) {
			if (_dims[d].type == H) {
				_init.push_back(h0[_dims[d].species1]);
			} else if (_dims[d].type == J) {
				swtch=true;
				if (j0.find(_dims[d].species1) != j0.end()) {
					if (j0[_dims[d].species1].find(_dims[d].species2) != j0[_dims[d].species1].end()) {
						_init.push_back(j0[_dims[d].species1][_dims[d].species2]);
						swtch=false;
					};
				};
				if (swtch) {
					_init.push_back(j0[_dims[d].species2][_dims[d].species1]);
				};
			};
		};

		/*
		std::cout << "Read in: ";
		for (int i=0; i<_init.size(); i++) {
			std::cout << _init[i] << " ";
		};
		std::cout << std::endl;
		*/

		f.close();
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
		for (int i=0; i<_var_nu_n_dim; i++) {
			r *= exp(-pow(_soln_traj[_var_time][_var_nu_dim_idxs[i]]-_dims[_var_nu_dim_idxs[i]].grid[_var_idxs[i]],2)/(2.*0.1))/sqrt(2.*M_PI*0.1);

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


