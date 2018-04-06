#include "../include/dynamic_boltzmann.hpp"
#include <iostream>
#include "../include/general.hpp"
#include <ctime>
#include <fstream>
#include <algorithm>

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Struct for specifying a dimension
	****************************************/

	Dim::Dim(std::string name, DimType type, std::string species, std::vector<std::string> basis_func_dims, double min, double max, int n, double init) : Dim(name,type,species, "", basis_func_dims, min, max, n, init) {};
	Dim::Dim(std::string name, DimType type, std::string species1, std::string species2, std::vector<std::string> basis_func_dims, double min, double max, int n, double init)
	{
		if ( (species1 == "") 
			|| 
			((type==H && species2 != "") || (type==W && species2 != "") || (type==J && species2 == ""))
			) {
			std::cerr << "ERROR! Dim specification is incorrect." << std::endl;
			exit(EXIT_FAILURE);
		};

		this->name = name;
		this->type = type;
		this->species1 = species1;
		this->species2 = species2;
		this->min = min;
		this->max = max;
		this->n = n;
		this->init = init;
		this->basis_func_dims = basis_func_dims;
	};

	/****************************************
	OptProblem
	****************************************/

	/********************
	Constructor
	********************/

	OptProblem::OptProblem(std::vector<Dim> dims, std::vector<std::string> species, double t_max, int n_t, int batch_size, int box_length, double dopt, int n_opt, int lattice_dim) : _latt(lattice_dim,box_length), _time("time",0.0,t_max,n_t)
	{
		// Set parameters
		if (DIAG_SETUP) { std::cout << "Copying params..." << std::flush; };
		_n_param = dims.size();
		_dopt = dopt;
		_n_opt = n_opt;
		_n_cd_steps = 1; // default
		_box_length = box_length;
		_n_batch = batch_size;
		_t_opt = 0;
		_n_t_soln = 0;
		_dir_io = "data/"; // default
		_fname_start_idx = 0; // default
		_write_bf_only_last = false; // default
		_hidden_layer_exists = false; // default
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Create the species and add to the lattice
		if (DIAG_SETUP) { std::cout << "Create species..." << std::flush; };
		for (auto s: species) {
			_species.push_back(Species(s));

			// Add to the lattice
			_latt.add_species(&(_species.back()));
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Create the interaction params
		if (DIAG_SETUP) { std::cout << "Create ixn params..." << std::flush; };
		for (auto d: dims) 
		{
			if (d.type==H) { 
				_ixn_params.push_back(IxnParamTraj(d.name,Hp,_find_species(d.species1),d.min,d.max,d.n,d.init,n_t));
			} else if (d.type==J) { 
				_ixn_params.push_back(IxnParamTraj(d.name,Jp,_find_species(d.species1),_find_species(d.species2),d.min,d.max,d.n,d.init,n_t));
			} else if (d.type==W) { 
				_ixn_params.push_back(IxnParamTraj(d.name,Wp,_find_species(d.species1),d.min,d.max,d.n,d.init,n_t));
			};
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Add the interaction params and time ptr to the species
		if (DIAG_SETUP) { std::cout << "Add ixn params to species..." << std::flush; };
		Species *sp1=nullptr, *sp2=nullptr;
		IxnParamTraj *ip_ptr=nullptr;
		for (auto d: dims) {
			ip_ptr = _find_ixn_param(d.name);
			if (d.type==H) {
				sp1 = _find_species(d.species1);
				sp1->set_h_ptr(ip_ptr);
			} else if (d.type==J) {
				sp1 = _find_species(d.species1);
				sp2 = _find_species(d.species2);
				sp1->add_j_ptr(sp2,ip_ptr);
				sp2->add_j_ptr(sp1,ip_ptr);
			} else if (d.type==W) {
				sp1 = _find_species(d.species1);
				sp1->set_w_ptr(ip_ptr);		
			};
		};
		for (auto itsp = _species.begin(); itsp!=_species.end(); itsp++) {
			itsp->set_opt_time_ptr(&_t_opt);
		};	
		// Ensure the J of the species are complete - if there is no ixn param that descripes the coupling, add a nullptr entry in the dictionary - later check if nullptr, then return 0
		// I think this is faster - otherwise there would be no reason to do it
		for (auto itsp1 = _species.begin(); itsp1!=_species.end(); itsp1++) {
			for (auto itsp2 = _species.begin(); itsp2!=_species.end(); itsp2++) {
				ip_ptr = _find_ixn_param_j_by_species(itsp1->name(), itsp2->name(), false);
				if (!ip_ptr) {
					// It's null; add to both
					itsp1->add_j_ptr(&(*itsp2),nullptr);
					itsp2->add_j_ptr(&(*itsp1),nullptr);
				};
			};
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Create the basis functions
		if (DIAG_SETUP) { std::cout << "Create basis funcs..." << std::flush; };
		std::vector<IxnParamTraj*> bf_ips;
		for (auto d: dims) {
			// Find the basis func dimensions
			bf_ips.clear();
			for (auto bfd: d.basis_func_dims) {
				bf_ips.push_back(_find_ixn_param(bfd));
			};
			// Make the basis function
			_bfs.push_back(BasisFunc("F_"+d.name,bf_ips));
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Add the basis functions to the interaction params
		if (DIAG_SETUP) { std::cout << "Add basis func to ixn params..." << std::flush; };
		BasisFunc* bf_ptr=nullptr;
		for (auto itp=_ixn_params.begin(); itp!=_ixn_params.end(); itp++) {
			// Find the basis func
			bf_ptr = _find_basis_func("F_"+itp->name());
			// Set
			itp->set_basis_func_ptr(bf_ptr);
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Create the variational terms
		if (DIAG_SETUP) { std::cout << "Create var terms..." << std::flush; };
		BasisFunc *num_bf_ptr = nullptr;
		IxnParamTraj *ixn_param_ptr = nullptr;
		// Go through numerators
		for (auto num: dims) {
			// Go through denominators
			for (auto denom: dims) {
				// Find the basis func dimensions
				bf_ips.clear();
				for (auto bfd: denom.basis_func_dims) {
					bf_ips.push_back(_find_ixn_param(bfd));
				};
				// Find the basis func
				bf_ptr = _find_basis_func("F_"+denom.name);
				// Find the interaction param
				ixn_param_ptr = _find_ixn_param(num.name);
				// Find the basis func corresponding to the numerator
				num_bf_ptr = _find_basis_func("F_"+num.name);
				// Create the var term
				_var_terms.push_back(VarTermTraj("var_"+num.name+"_wrt_F_"+denom.name, ixn_param_ptr, bf_ptr, bf_ips, num_bf_ptr, num.basis_func_dims.size(), n_t));
			};
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Add the pointers to variational terms needed to update this var term
		if (DIAG_SETUP) { std::cout << "Add ptrs to var term..." << std::flush; };
		VarTermTraj* vt_ptr=nullptr;
		// Go through numerators
		for (auto num: dims) {
			// Go through denoms
			for (auto denom=_bfs.begin(); denom!=_bfs.end(); denom++) {
				// Find the ixn params that are arguments to the num's basis func
				bf_ips.clear();
				for (auto bfd: num.basis_func_dims) {
					bf_ips.push_back(_find_ixn_param(bfd));
				};
				// Find the variational term
				vt_ptr = _find_var_term("var_"+num.name+"_wrt_"+denom->name());
				// Find the variational terms needed to update this one
				for (auto ip_ptr: bf_ips) {
					vt_ptr->add_update_ptr(_find_var_term("var_"+ip_ptr->name()+"_wrt_"+denom->name()));
				};
			};
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Add pointers to the basis functions needed to update them
		if (DIAG_SETUP) { std::cout << "Add ptrs to basis func..." << std::flush; };
		for (auto itbf=_bfs.begin(); itbf!=_bfs.end(); itbf++) {
			for (auto itp=_ixn_params.begin(); itp!=_ixn_params.end(); itp++) {
				// Find the variational term
				vt_ptr = _find_var_term("var_"+itp->name()+"_wrt_"+itbf->name());
				// Add
				itbf->add_update_ptrs(&*itp,vt_ptr);
			};
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };
	};

	OptProblem::OptProblem(const OptProblem& other) : _time(other._time)
	{
		_copy(other);
	};

	OptProblem::OptProblem(OptProblem&& other) : _time(other._time)
	{
		_copy(other);
		other._reset();
	};

	OptProblem& OptProblem::operator=(const OptProblem& other)
	{
		if (this != &other)
		{			
			_clean_up();
			_copy(other);
		};
		return *this;		
	};

	OptProblem& OptProblem::operator=(OptProblem&& other)
	{
		if (this != &other)
		{
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;		
	};

	OptProblem::~OptProblem() {
		_clean_up();
	};

	/********************
	Helpers for constructors
	********************/

	void OptProblem::_clean_up()
	{
		// Nothing...
	};
	void OptProblem::_reset()
	{
		_dir_io = "";
		_n_param = 0;
		_ixn_params.clear();
		_bfs.clear();
		_var_terms.clear();
		_hidden_layer_exists = false;
		_hidden_units.clear();
		_species.clear();
		_fnames.clear();
		_n_t_soln = 0;
		_t_opt = 0;
		_n_batch = 0;
		_n_cd_steps = 1;
		_box_length = 0;
		_dopt = 0;
		_fname_start_idx = 0;
		_write_bf_only_last = false;
	};

	void OptProblem::_copy(const OptProblem& other)
	{
		_dir_io = other._dir_io;
		_n_param = other._n_param;
		_ixn_params = other._ixn_params;
		_bfs = other._bfs;
		_var_terms = other._var_terms;
		_hidden_layer_exists = other._hidden_layer_exists;
		_hidden_units = other._hidden_units;
		_time = other._time;
		_species = other._species;
		_fnames = other._fnames;
		_n_t_soln = other._n_t_soln;
		_t_opt = other._t_opt;
		_n_batch = other._n_batch;
		_n_cd_steps = other._n_cd_steps;
		_box_length = other._box_length;
		_latt = other._latt;
		_dopt = other._dopt;
		_n_opt = other._n_opt;
		_fname_start_idx = other._fname_start_idx;
		_write_bf_only_last = other._write_bf_only_last;
	};

	/********************
	Set properties
	********************/

	void OptProblem::add_hidden_unit(std::vector<std::vector<int>> lattice_idxs, std::string species) 
	{
		// Find sites indicated by connections
		std::vector<Site*> conns;
		for (auto c: lattice_idxs) {
			if (c.size() != _latt.dim()) {
				std::cerr << "ERROR: lattice dim does not match hidden layer dim." << std::endl;
				exit(EXIT_FAILURE);
			};
			if (c.size() == 1) {
				conns.push_back(_latt.get_site(c[0]));
			} else if (c.size() == 2) {
				conns.push_back(_latt.get_site(c[0],c[1]));
			} else if (c.size() == 3) {
				conns.push_back(_latt.get_site(c[0],c[1],c[2]));
			};
		};

		// Add
		_add_hidden_unit(conns,species);
	};
	void OptProblem::add_hidden_unit(std::vector<int> lattice_idxs, std::string species) 
	{
		// Find sites indicated by connections
		std::vector<Site*> conns;
		for (auto c: lattice_idxs) {
			conns.push_back(_latt.get_site(c));
		};

		// Add
		_add_hidden_unit(conns,species);
	};
	void OptProblem::_add_hidden_unit(std::vector<Site*> conns, std::string species)
	{
		// Flag that a hidden layer exists
		_hidden_layer_exists = true;
		// Also tell the lattice - important for the annealer
		_latt.set_hidden_layer_exists();

		// Find the species
		Species *sp = _find_species(species);

		// Make hidden unit
		_hidden_units.push_back(HiddenUnit(conns,sp));

		// Go through lattice sites in this connection
		// Indicate that for this species, they are linked to this hidden unit
		for (auto s: conns) {
			s->hidden_conns[sp].push_back(&_hidden_units.back());
		};

		// Tell the appropriate interaction parameter that these these sites are connected to this hidden unit
		IxnParamTraj *ip = _find_ixn_param_w_by_species(species);
		for (auto sptr: conns) {
			ip->add_visible_hidden_connection(sptr,&_hidden_units.back());
		};
	};

	void OptProblem::set_dir_io(std::string dir) {
		_dir_io = dir;
	};

	void OptProblem::set_fname_start_idx(int idx) {
		_fname_start_idx = idx;
	};

	void OptProblem::add_fname(std::string f) {
		_fnames.push_back(f);
	};

	void OptProblem::set_n_cd_steps(int n_steps) {
		_n_cd_steps = n_steps;
	};

	/********************
	Validate setup
	********************/

	void OptProblem::validate_setup() const {
		std::cout << "------------------------" << std::endl;
		for (auto it: _species) {
			it.validate_setup();
		};
		for (auto it: _ixn_params) {
			it.validate_setup();
		};
		for (auto it: _var_terms) {
			it.validate_setup();
		};
		for (auto it: _bfs) {
			it.validate_setup();
		};
		std::cout << "------------------------" << std::endl;
	};

	/********************
	Solve interaction parameter traj
	********************/

	void OptProblem::solve_ixn_param_traj() {
		// Number of time points before the solution has gone out of bounds
		_n_t_soln = 1;

		// Have we gone out of the domain?
		bool in_domain;

		// Go through all times
		for (int it=1; it<_time.n(); it++) {
			// Go through all ixn params
			for (auto itp=_ixn_params.begin(); itp!=_ixn_params.end(); itp++) {
				// Solve
				in_domain = itp->calculate_at_time(it,_time.delta());
				if (!in_domain) {
					// Stop
					std::cout << "WARNING: out of bounds: " << itp->get_at_time(it) << " is out of the grid at time: " << _n_t_soln+1 << std::endl;
					itp->print_grid_range();
					return;
				};
			};
			// Increment time
			_n_t_soln++;
		};
	};

	/********************
	Solve variational term traj
	********************/

	void OptProblem::solve_var_traj() {
		// Go through all times
		for (int it=1; it<_n_t_soln; it++) {
			// Go through all var terms
			for (auto itv=_var_terms.begin(); itv!=_var_terms.end(); itv++) {
				//std::cout << "calc " << itv->name() << " at time " << it << std::endl;
				itv->calculate_at_time(it,_time.delta());
				//std::cout << "ok" << std::endl;
			};
		};
	};

	/********************
	Solve --- Main optimization loop
	********************/

	void OptProblem::solve(bool verbose, bool same_lattice)
	{
		// Write the grids
		write_bf_grids();
		write_t_grid();

		// Iterate over optimization steps
		for (int i_opt=0; i_opt<_n_opt; i_opt++)
		{
			std::cout << "Opt step " << i_opt << " / " << _n_opt-1 << std::endl;

			// Write the basis funcs
			if (!_write_bf_only_last) {
				write_bfs(_dir_io+"F/",i_opt);
			};

			/*****
			Step 1 - Solve the current trajectory
			*****/

			if (DIAG_SOLVE) { std::cout << "Solving ixn param" << std::endl; };

			solve_ixn_param_traj();

			// Write
			write_ixn_params(_dir_io+"ixn_params/",i_opt);

			if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

			/*****
			Step 2 - Solve the variational problem traj
			*****/

			if (DIAG_SOLVE) { std::cout << "Solving var term" << std::endl; };

			solve_var_traj();

			// Write
			//write_var_terms(_dir_io+"var_terms/",i_opt);
			
			if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

			/*****
			Step 3 - Pick a random batch
			*****/

			if (DIAG_SOLVE) { std::cout << "Random batch" << std::endl; };

			std::vector<std::string> fnames, fnames_possible=_fnames;
			std::vector<std::string>::iterator itf;
			if (same_lattice == false) {
				for (int i_batch=0; i_batch<_n_batch; i_batch++) {
					itf = fnames_possible.begin();
					std::advance(itf,randI(0,fnames_possible.size()-1));
					fnames.push_back(*itf);
					fnames_possible.erase(itf); // don't choose again
				};
			} else {
				// Pick a rand
				itf = fnames_possible.begin();
				std::advance(itf,randI(0,fnames_possible.size()-1));
				fnames.push_back(*itf);
				for (int i_batch=1; i_batch<_n_batch; i_batch++) {
					// Only this one
					fnames.push_back(fnames.back());
				};
			};

			if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

			/*****
			Step 4 - reset the moments at all times
			*****/

			if (DIAG_SOLVE) { std::cout << "Reset moments" << std::endl; };

			for (auto itp = _ixn_params.begin(); itp != _ixn_params.end(); itp++) {
				itp->moments_reset();
			};

			if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

			/*****
			Step 5 - Go through all times
			*****/

			if (DIAG_SOLVE) { std::cout << "Looping over times" << std::endl; };

			// Use the class variable _t_opt to iterate
			for (_t_opt=0; _t_opt < _n_t_soln; _t_opt++)
			{
				if (verbose) {
					std::cout << "time: " << _t_opt << std::flush;
				};

				/*****
				Step 5.1 - loop over all samples in the batch
				*****/

				if (DIAG_SOLVE) { std::cout << "   Looping over batch" << std::endl; };

				for (int i_batch=0; i_batch<_n_batch; i_batch++) 
				{
					if (verbose) {
						std::cout << "." << std::flush;
					};

					/*****
					Step 5.1.1 - Read in batch at this timestep
					*****/

					if (DIAG_SOLVE) { std::cout << "      Read in batch" << std::endl; };

					_latt.read_from_file(fnames[i_batch] + pad_str(_fname_start_idx+_t_opt,4) + ".txt");

					/*****
					Step 5.1.2 - Activate the hidden units if needed
					*****/

					if (DIAG_SOLVE) { std::cout << "      Activating hidden units" << std::endl; };

					if (_hidden_layer_exists) {
						for (auto ithu = _hidden_units.begin(); ithu != _hidden_units.end(); ithu++) {
							// When using real data, always use binary states
							ithu->activate(false);
						};
					};

					/*****
					Step 5.1.3 - Record the awake moments
					*****/

					if (DIAG_SOLVE) { std::cout << "      Record awake moments" << std::endl; };

					for (auto itp = _ixn_params.begin(); itp != _ixn_params.end(); itp++) {
						itp->moments_retrieve_at_time(IxnParamTraj::AWAKE,_t_opt,_n_batch);
					};

					/*****
					Step 5.1.4 - CD steps
					*****/

					for (int cd_step=0; cd_step<_n_cd_steps; cd_step++)
					{
						// Sample

						if (DIAG_SOLVE) { std::cout << "      Anneal" << std::endl; };

						_latt.sample();

						// Activate the hidden units if needed

						if (DIAG_SOLVE) { std::cout << "      Activating hidden units" << std::endl; };

						if (_hidden_layer_exists) {
							for (auto ithu = _hidden_units.begin(); ithu != _hidden_units.end(); ithu++) {
								// When using reconstructions, always use raw probabilities
								ithu->activate(false);
							};
						};
					};

					/*****
					Step 5.1.5 - Record the asleep moments
					*****/

					if (DIAG_SOLVE) { std::cout << "      Record asleep moments" << std::endl; };

					for (auto itp = _ixn_params.begin(); itp != _ixn_params.end(); itp++) {
						itp->moments_retrieve_at_time(IxnParamTraj::ASLEEP,_t_opt,_n_batch);
					};
				};

				if (verbose) {
					std::cout << std::endl;
				};

				if (DIAG_SOLVE) { std::cout << "   OK" << std::endl; };
			};

			if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

			// Write the moments
			write_moments(_dir_io+"moments/",i_opt);

			/*****
			Step 6 - Update the basis funcs
			*****/

			for (auto itbf=_bfs.begin(); itbf!=_bfs.end(); itbf++) {
				itbf->update(_n_t_soln, _time.delta(), _dopt);
			};
		};

		// Write the basis funcs one last time
		write_bfs(_dir_io+"F/",_n_opt);
	};



	/********************
	Solve over varying initial conditions
	********************/

	void OptProblem::solve_varying_ic(bool verbose)
	{
		// Write the grids
		write_bf_grids();
		write_t_grid();

		// Declare
		std::vector<std::string> fnames, fnames_possible;
		std::vector<std::string>::iterator itf;
		std::vector<int> batch_idxs;

		// Iterate over optimization steps
		for (int i_opt=0; i_opt<_n_opt; i_opt++)
		{
			std::cout << "Opt step " << i_opt << " / " << _n_opt-1 << std::endl;

			// Write the basis funcs
			write_bfs(_dir_io+"F/",i_opt);

			/*****
			Step 1 - Pick a random batch
			*****/

			if (DIAG_SOLVE) { std::cout << "Random batch" << std::endl; };

			std::cout << "Random batch" << std::endl;
			fnames.clear();
			fnames_possible=_fnames;
			batch_idxs.clear();
			for (int i_batch=0; i_batch<_n_batch; i_batch++) {
				itf = fnames_possible.begin();
				std::advance(itf,randI(0,fnames_possible.size()-1));
				fnames.push_back(*itf);
				batch_idxs.push_back(find(_fnames.begin(), _fnames.end(), fnames.back()) - _fnames.begin());
				fnames_possible.erase(itf);
			};
			std::cout << "OK" << std::endl;

			if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

			/*****
			Step 2 - Iterate over batch
			*****/

			if (DIAG_SOLVE) { std::cout << "   Looping over batch" << std::endl; };

			for (int i_batch=0; i_batch<_n_batch; i_batch++) 
			{
				if (verbose) {
					std::cout << "sample: " << i_batch << " / " << _n_batch << std::endl;
				};

				/*****
				Step 2.1 - Read the IC
				*****/

				read_init_cond(fnames[i_batch]+"../");

				/*****
				Step 2.2 - Solve the current trajectory
				*****/

				if (DIAG_SOLVE) { std::cout << "Solving ixn param" << std::endl; };

				solve_ixn_param_traj();

				// Write
				write_ixn_params(_dir_io+"ixn_params/",i_opt,batch_idxs[i_batch]);

				if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

				/*****
				Step 2.3 - Solve the variational problem traj
				*****/

				if (DIAG_SOLVE) { std::cout << "Solving var term" << std::endl; };

				solve_var_traj();

				// Write
				//write_var_terms(_dir_io+"var_terms/",i_opt);
				
				if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

				/*****
				Step 2.4 - reset the moments at all times
				*****/

				if (DIAG_SOLVE) { std::cout << "Reset moments" << std::endl; };

				for (auto itp = _ixn_params.begin(); itp != _ixn_params.end(); itp++) {
					itp->moments_reset();
				};

				if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

				/*****
				Step 2.5 - Go through all times
				*****/

				if (DIAG_SOLVE) { std::cout << "Looping over times" << std::endl; };

				// Use the class variable _t_opt to iterate
				for (_t_opt=0; _t_opt < _n_t_soln; _t_opt++)
				{
					if (verbose) {
						std::cout << "." << std::flush;
					};

					/*****
					Step 2.5.1 - Read in sample at this timestep
					*****/

					if (DIAG_SOLVE) { std::cout << "      Read in batch" << std::endl; };

					_latt.read_from_file(fnames[i_batch] + pad_str(_fname_start_idx+_t_opt,4) + ".txt");

					/*****
					Step 2.5.2 - Record the awake moments
					*****/

					if (DIAG_SOLVE) { std::cout << "      Record awake moments" << std::endl; };

					for (auto itp = _ixn_params.begin(); itp != _ixn_params.end(); itp++) {
						itp->moments_retrieve_at_time(IxnParamTraj::AWAKE,_t_opt);
					};

					/*****
					Step 2.5.3 - CD steps - alternate sampling and activating
					*****/

					for (int cd_step=0; cd_step<_n_cd_steps; cd_step++) {

						// Anneal

						if (DIAG_SOLVE) { std::cout << "      Anneal" << std::endl; };

						_latt.sample();

						// Record the asleep moments

						if (DIAG_SOLVE) { std::cout << "      Record asleep moments" << std::endl; };

						for (auto itp = _ixn_params.begin(); itp != _ixn_params.end(); itp++) {
							itp->moments_retrieve_at_time(IxnParamTraj::ASLEEP,_t_opt);
						};

					};
				};

				if (verbose) {
					std::cout << std::endl;
				};

				if (DIAG_SOLVE) { std::cout << "   OK" << std::endl; };

				/*****
				Step 2.6 - Write the moments
				*****/

				write_moments(_dir_io+"moments/",i_opt,batch_idxs[i_batch]);

				/*****
				Step 2.7 - Gather the update (but dont commit)
				*****/

				for (auto itbf=_bfs.begin(); itbf!=_bfs.end(); itbf++) {
					itbf->update_gather(_n_t_soln, _time.delta(), _dopt);
				};

			};

			/*****
			Step 3 - Commit the updates to the basis funcs
			*****/

			for (auto itbf=_bfs.begin(); itbf!=_bfs.end(); itbf++) {
				itbf->update_committ_gathered();
			};

		};

		// Write the basis funcs one last time
		write_bfs(_dir_io+"F/",_n_opt);
	};


	/********************
	Read some initial conditions
	********************/

	void OptProblem::read_init_cond(std::string dir) {

		std::ifstream f;
		f.open(dir+"init.txt", std::ios::in);
		char frag[100]; // fragments of the line
		std::string sname="",sval="";
		int i_frag=0;
		IxnParamTraj *ip;

		if (f.is_open()) { // make sure we found it
			while (!f.eof()) {
				f >> frag;
				if (i_frag==0) {
					sname += frag; i_frag++;
				} else if (i_frag==1) {
					sval += frag;

					// Find the ixn param with this name
					ip = _find_ixn_param(sname);
					if (ip) {
						ip->set_init_cond(std::stod(sval));
					} else {
						std::cerr << "ERROR: Ixn param with name " << sname << " not found while reading IC." << std::endl;
						exit(EXIT_FAILURE); 
					};

					// Reset
					sname=""; sval=""; i_frag=0;
				};
			};
		};

		f.close();
	};

	/********************
	Writing functions
	********************/

	void OptProblem::write_bf_grids() const {
		for (auto it=_bfs.begin(); it!=_bfs.end(); it++) {
			it->write_grid(_dir_io+"grid_"+it->name()+".txt");
		};
	};
	void OptProblem::write_t_grid() const {
		_time.write_grid(_dir_io+"grid_time.txt");
	};

	void OptProblem::write_ixn_params(std::string dir, int idx) const {
		for (auto it = _ixn_params.begin(); it!=_ixn_params.end(); it++) {
			it->write_vals(dir,idx,_n_t_soln);
		};
	};
	void OptProblem::write_ixn_params(std::string dir, int idx1, int idx2) const {
		for (auto it = _ixn_params.begin(); it!=_ixn_params.end(); it++) {
			it->write_vals(dir,idx1,idx2,_n_t_soln);
		};
	};
	void OptProblem::write_bfs(std::string dir, int idx) const {
		for (auto it=_bfs.begin(); it!=_bfs.end(); it++) {
			it->write_vals(dir, idx);
		};
	};
	void OptProblem::write_var_terms(std::string dir, int idx) const {
		for (auto it=_var_terms.begin(); it!=_var_terms.end(); it++) {
			it->write_vals(dir,idx);
		};
	};
	void OptProblem::write_moments(std::string dir, int idx) const {
		for (auto it = _ixn_params.begin(); it!=_ixn_params.end(); it++) {
			it->write_moments(dir,idx,_n_t_soln);
		};
	};
	void OptProblem::write_moments(std::string dir, int idx1, int idx2) const {
		for (auto it = _ixn_params.begin(); it!=_ixn_params.end(); it++) {
			it->write_moments(dir,idx1,idx2,_n_t_soln);
		};
	};

	void OptProblem::set_flag_write_bf_only_final()
	{
		_write_bf_only_last = true;
	};

	/********************
	Read
	********************/

	void OptProblem::read_bf(std::string bf_name, std::string fname) 
	{
		// Find the basis func
		BasisFunc* bf = _find_basis_func(bf_name);
		if (!bf) {
			std::cerr << "ERROR: Could not find basis func to read into." << std::endl;
			exit(EXIT_FAILURE);
		};
		bf->read_vals(fname);
	};

	/****************************************
	OptProblem - PRIVATE
	****************************************/

	/********************
	Search functions
	********************/

	Species* OptProblem::_find_species(std::string name, bool enforce_success) {
		for (auto it=_species.begin(); it!=_species.end(); it++) {
			if (it->name() == name) {
				return &*it;
			};
		};
		if (enforce_success) {
			std::cerr << "ERROR: could not find species: " << name << std::endl;
			exit(EXIT_FAILURE);
		} else {
			return nullptr;
		};
	};
	IxnParamTraj* OptProblem::_find_ixn_param(std::string name, bool enforce_success) {
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			if (it->name() == name) {
				return &*it;
			};
		};
		if (enforce_success) {
			std::cerr << "ERROR: could not find ixn param: " << name << std::endl;
			exit(EXIT_FAILURE);
		} else {
			return nullptr;
		};
	};
	IxnParamTraj* OptProblem::_find_ixn_param_j_by_species(std::string species_name_1, std::string species_name_2, bool enforce_success) {
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			if (it->is_j_with_species(species_name_1,species_name_2)) {
				return &*it;
			};
		};
		if (enforce_success) {
			std::cerr << "ERROR: could not find J ixn param for species: " << species_name_1 << " " << species_name_2 << std::endl;
			exit(EXIT_FAILURE);
		} else {
			return nullptr;
		};
	};
	IxnParamTraj* OptProblem::_find_ixn_param_w_by_species(std::string species_name, bool enforce_success) {
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			if (it->is_w_with_species(species_name)) {
				return &*it;
			};
		};
		if (enforce_success) {
			std::cerr << "ERROR: could not find visible-to-hidden ixn param for species: " << species_name << std::endl;
			exit(EXIT_FAILURE);
		} else {
			return nullptr;
		};
	};
	BasisFunc* OptProblem::_find_basis_func(std::string name, bool enforce_success) {
		for (auto it=_bfs.begin(); it!=_bfs.end(); it++) {
			if (it->name() == name) {
				return &*it;
			};
		};
		if (enforce_success) {
			std::cerr << "ERROR: could not find basis func: " << name << std::endl;
			exit(EXIT_FAILURE);
		} else {
			return nullptr;
		};
	};
	VarTermTraj* OptProblem::_find_var_term(std::string name, bool enforce_success) {
		for (auto it=_var_terms.begin(); it!=_var_terms.end(); it++) {
			if (it->name() == name) {
				return &*it;
			};
		};
		if (enforce_success) {
			std::cerr << "ERROR: could not find var term: " << name << std::endl;
			exit(EXIT_FAILURE);
		} else {
			return nullptr;
		};
	};
};