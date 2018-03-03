#include "dynamic_boltzmann.hpp"
#include <iostream>
#include "general.hpp"
#include <ctime>

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
		if ((type==H && species2 != "") || (type==J && species2 == "")) {
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

	OptProblem::OptProblem(std::vector<Dim> dims, std::vector<std::string> species, double t_max, int n_t, int batch_size, int n_annealing, int box_length, double dopt, int n_opt) : _latt(box_length), _time("time",0.0,t_max,n_t)
	{
		// Set parameters
		if (DIAG_SETUP) { std::cout << "Copying params..." << std::flush; };
		_n_param = dims.size();
		_dopt = dopt;
		_n_opt = n_opt;
		_n_annealing = n_annealing;
		_box_length = box_length;
		_n_batch = batch_size;
		_t_opt = 0;
		_n_t_soln = 0;
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
		for (auto d: dims) {
			if (d.type==H) { 
				_ixn_params.push_back(IxnParam(d.name,Hp,_find_species(d.species1),d.min,d.max,d.n,d.init,n_t));
			} else { 
				_ixn_params.push_back(IxnParam(d.name,Jp,_find_species(d.species1),_find_species(d.species2),d.min,d.max,d.n,d.init,n_t));
			};
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Add the interaction params and time ptr to the species
		if (DIAG_SETUP) { std::cout << "Add ixn params to species..." << std::flush; };
		Species *sp1=nullptr, *sp2=nullptr;
		IxnParam *ip_ptr=nullptr;
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
			};
		};
		for (auto itsp = _species.begin(); itsp!=_species.end(); itsp++) {
			itsp->set_opt_time_ptr(&_t_opt);
		};	
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Create the basis functions
		if (DIAG_SETUP) { std::cout << "Create basis funcs..." << std::flush; };
		std::vector<IxnParam*> bf_ips;
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
		IxnParam *ixn_param_ptr = nullptr;
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
				_var_terms.push_back(VarTerm("var_"+num.name+"_wrt_F_"+denom.name, ixn_param_ptr, bf_ptr, bf_ips, num_bf_ptr, num.basis_func_dims.size(), n_t));
			};
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Add the pointers to variational terms needed to update this var term
		if (DIAG_SETUP) { std::cout << "Add ptrs to var term..." << std::flush; };
		VarTerm* vt_ptr=nullptr;
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
	OptProblem::~OptProblem() {
		// Nothing...
	};

	/********************
	Set properties
	********************/

	void OptProblem::add_fname(std::string f) {
		_fnames.push_back(f);
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

	void OptProblem::solve(bool verbose)
	{
		// Write the grids
		write_bf_grids();
		write_t_grid();

		// Iterate over optimization steps
		for (int i_opt=0; i_opt<_n_opt; i_opt++)
		{
			std::cout << "Opt step " << i_opt << " / " << _n_opt-1 << std::endl;

			// Write the basis funcs
			write_bfs("data/F/",i_opt);

			/*****
			Step 1 - Solve the current trajectory
			*****/

			if (DIAG_SOLVE) { std::cout << "Solving ixn param" << std::endl; };

			solve_ixn_param_traj();

			// Write
			write_ixn_params("data/ixn_params/",i_opt);

			if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

			/*****
			Step 2 - Solve the variational problem traj
			*****/

			if (DIAG_SOLVE) { std::cout << "Solving var term" << std::endl; };

			solve_var_traj();

			// Write
			//write_var_terms("data/var_terms/",i_opt);
			
			if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

			/*****
			Step 3 - Pick a random batch
			*****/

			if (DIAG_SOLVE) { std::cout << "Random batch" << std::endl; };

			std::vector<std::string> fnames, fnames_possible=_fnames;
			std::vector<std::string>::iterator itf;
			for (int i_batch=0; i_batch<_n_batch; i_batch++) {
				itf = fnames_possible.begin();
				std::advance(itf,randI(0,fnames_possible.size()-1));
				fnames.push_back(*itf);
				// fnames_possible.erase(itf);
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

					_latt.read_from_file(fnames[i_batch] + pad_str(_t_opt,4) + ".txt");

					/*****
					Step 5.1.2 - Record the awake moments
					*****/

					if (DIAG_SOLVE) { std::cout << "      Record awake moments" << std::endl; };

					for (auto itp = _ixn_params.begin(); itp != _ixn_params.end(); itp++) {
						itp->moments_retrieve_at_time(IxnParam::AWAKE,_t_opt,_n_batch);
					};

					/*****
					Step 5.1.3 - anneal
					*****/

					if (DIAG_SOLVE) { std::cout << "      Anneal" << std::endl; };

					_latt.anneal(_n_annealing);

					/*****
					Step 5.1.4 - Record the asleep moments
					*****/

					if (DIAG_SOLVE) { std::cout << "      Record asleep moments" << std::endl; };

					for (auto itp = _ixn_params.begin(); itp != _ixn_params.end(); itp++) {
						itp->moments_retrieve_at_time(IxnParam::ASLEEP,_t_opt,_n_batch);
					};
				};

				if (verbose) {
					std::cout << std::endl;
				};

				if (DIAG_SOLVE) { std::cout << "   OK" << std::endl; };
			};

			if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

			// Write the moments
			write_moments("data/moments/",i_opt);

			/*****
			Step 6 - Update the basis funcs
			*****/

			for (auto itbf=_bfs.begin(); itbf!=_bfs.end(); itbf++) {
				itbf->update(_n_t_soln, _time.delta(), _dopt);
			};
		};

		// Write the basis funcs one last time
		write_bfs("data/F/",_n_opt);
	};

	/********************
	Writing functions
	********************/

	void OptProblem::write_bf_grids() const {
		for (auto it=_bfs.begin(); it!=_bfs.end(); it++) {
			it->write_grid("data/grid_"+it->name()+".txt");
		};
	};
	void OptProblem::write_t_grid() const {
		_time.write_grid("data/grid_time.txt");
	};

	void OptProblem::write_ixn_params(std::string dir, int idx) const {
		for (auto it = _ixn_params.begin(); it!=_ixn_params.end(); it++) {
			it->write_vals(dir,idx,_n_t_soln);
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

	/****************************************
	OptProblem - PRIVATE
	****************************************/

	/********************
	Search functions
	********************/

	Species* OptProblem::_find_species(std::string name) {
		for (auto it=_species.begin(); it!=_species.end(); it++) {
			if (it->name() == name) {
				return &*it;
			};
		};
		std::cerr << "ERROR: could not find species: " << name << std::endl;
		exit(EXIT_FAILURE);
	};
	IxnParam* OptProblem::_find_ixn_param(std::string name) {
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			if (it->name() == name) {
				return &*it;
			};
		};
		std::cerr << "ERROR: could not find ixn param: " << name << std::endl;
		exit(EXIT_FAILURE);
	};
	BasisFunc* OptProblem::_find_basis_func(std::string name) {
		for (auto it=_bfs.begin(); it!=_bfs.end(); it++) {
			if (it->name() == name) {
				return &*it;
			};
		};
		std::cerr << "ERROR: could not find basis func: " << name << std::endl;
		exit(EXIT_FAILURE);
	};
	VarTerm* OptProblem::_find_var_term(std::string name) {
		for (auto it=_var_terms.begin(); it!=_var_terms.end(); it++) {
			if (it->name() == name) {
				return &*it;
			};
		};
		std::cerr << "ERROR: could not find var term: " << name << std::endl;
		exit(EXIT_FAILURE);
	};

};