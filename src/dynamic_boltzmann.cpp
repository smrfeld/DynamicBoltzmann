#include "../include/dynamic_boltzmann.hpp"
#include <iostream>
#include "../include/general.hpp"
#include <ctime>
#include <fstream>
#include <algorithm>
#include <set>
#include "var_term_traj.hpp"
#include <stdlib.h>

#define DIAG_SETUP 0
#define DIAG_SOLVE 0

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
		if (type == B || type == H || type == W) {
			if (species1 == "" || species2 != "") {
				std::cerr << "ERROR! Dim specification is incorrect for H, B or W." << std::endl;
				exit(EXIT_FAILURE);
			};
		} else if (type == J) {
			if (species1 == "" || species2 == "") {
				std::cerr << "ERROR! Dim specification is incorrect for J." << std::endl;
				exit(EXIT_FAILURE);
			};
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
	OptProblem - IMPLEMENTATION
	****************************************/

	class OptProblem::Impl {

	private:

		/********************
		Parameters
		********************/

		// Number of dimensions
		int _n_param;

		// List of interaction parameters
		std::list<IxnParamTraj> _ixn_params;

		// List of basis funcs
		std::list<BasisFunc> _bfs;

		// List of variational terms
		std::list<VarTermTraj> _var_terms;

		// List of hidden units, and flag if they exist
		bool _hidden_layer_exists;
		std::list<HiddenUnit> _hidden_units;

		// Time dimension
		Grid _time;

		// Species present
		std::list<Species> _species;

		// Number of steps in this nu solution
		int _n_t_soln;

		// The current time in the optimization step
		int _t_opt;

		// Lattice size
		int _box_length;

		// Lattice to hold the current sample of the batch
		Lattice _latt;

		// Add a hidden unit
		void _add_hidden_unit(std::vector<Site*> conns, std::string species);

		// Search functions
		Species* _find_species(std::string name, bool enforce_success=true);
		IxnParamTraj* _find_ixn_param_by_name(std::string name, bool enforce_success=true);
		IxnParamTraj* _find_ixn_param_b_by_species(std::string species_name, bool enforce_success=true);
		IxnParamTraj* _find_ixn_param_j_by_species(std::string species_name_1, std::string species_name_2, bool enforce_success=true);
		IxnParamTraj* _find_ixn_param_w_by_species(std::string species_name, bool enforce_success=true);
		BasisFunc* _find_basis_func(std::string name, bool enforce_success=true);
		VarTermTraj* _find_var_term(std::string name, bool enforce_success=true);

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const Impl& other);

	public:

		/********************
		Constructor
		********************/

		Impl(std::vector<Dim> dims, double t_max, int n_t, int box_length, int lattice_dim);
		Impl(Impl&& other);
	    Impl& operator=(Impl&& other);
		~Impl();

		/********************
		Set properties	
		********************/

		// Any dim
		void add_hidden_unit(std::vector<std::vector<int>> lattice_idxs, std::string species);
		// 1D specific
		void add_hidden_unit(std::vector<int> lattice_idxs, std::string species);

		/********************
		Validate setup by printing
		********************/

		void validate_setup() const;

		/********************
		Solve interaction parameter traj
		********************/

		void solve_ixn_param_traj();

		/********************
		Solve for variational trajectory
		********************/

		void solve_var_traj();

		/********************
		Solve
		********************/

		void solve(std::vector<std::string> fnames, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolve options);
		void solve_varying_ic(std::vector<std::pair<std::string,int>> fnames, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolve options);

		/********************
		Read basis func
		********************/

		void read_basis_func(std::string bf_name, std::string fname);

		/********************
		Read some initial conditions
		********************/

		void read_init_cond(std::string dir);

		/********************
		Write
		********************/

		void write_bf_grids(std::string dir) const;
		void write_t_grid(std::string dir) const;

		void write_ixn_params(std::string dir, int idx) const;
		void write_ixn_params(std::string dir, int idx1, int idx2) const;
		void write_bfs(std::string dir, int idx) const;
		void write_var_terms(std::string dir, int idx) const;
		void write_moments(std::string dir, int idx) const;
		void write_moments(std::string dir, int idx1, int idx2) const;
	};

	/****************************************
	OptProblem - IMPLEMENTATION DEFINITIONS
	****************************************/

	/********************
	Constructor
	********************/

	OptProblem::Impl::Impl(std::vector<Dim> dims, double t_max, int n_t, int box_length, int lattice_dim) : _latt(lattice_dim,box_length), _time("time",0.0,t_max,n_t)
	{
		// Set parameters
		if (DIAG_SETUP) { std::cout << "Copying params..." << std::flush; };
		_n_param = dims.size();
		_box_length = box_length;
		_t_opt = 0;
		_n_t_soln = 0;
		_hidden_layer_exists = false;
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Set of species that exist
		if (DIAG_SETUP) { std::cout << "Creating set of species..." << std::flush; };
		std::set<std::string> species;
		for (auto d: dims) {
			if (d.species1 != "") { species.insert(d.species1); };
			if (d.species2 != "") { species.insert(d.species2); };
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Tell the lattice about what dims exist
		if (DIAG_SETUP) { std::cout << "Telling lattice what dims exist..." << std::flush; };
		for (auto d: dims) {
			if (d.type==H) {
				_latt.set_exists_h(true);
			} else if (d.type==J) {
				_latt.set_exists_j(true);
			} else if (d.type==W) {
				_hidden_layer_exists = true;
				_latt.set_exists_w(true);
			} else if (d.type==B) {
				_hidden_layer_exists = true;
				// No need to tell lattice, since it only affects hidden unit activation and not sampling
			};
		};
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
			} else if (d.type==B) {
				_ixn_params.push_back(IxnParamTraj(d.name,Bp,_find_species(d.species1),d.min,d.max,d.n,d.init,n_t));
			};
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Add the interaction params and time ptr to the species
		if (DIAG_SETUP) { std::cout << "Add ixn params to species..." << std::flush; };
		Species *sp1=nullptr, *sp2=nullptr;
		IxnParamTraj *ip_ptr=nullptr;
		for (auto d: dims) {
			ip_ptr = _find_ixn_param_by_name(d.name);
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
			// No need to tell species about biases
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

		// Initialize counting structures on species
		if (DIAG_SETUP) { std::cout << "Init counting on species..." << std::flush; };
		std::vector<Species*> sp_vec;
		for (auto it = _species.begin(); it != _species.end(); it++) {
			sp_vec.push_back(&(*it));
		};
		for (auto it = _species.begin(); it != _species.end(); it++) {
			it->count_nn_for_species(sp_vec);
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Create the basis functions
		if (DIAG_SETUP) { std::cout << "Create basis funcs..." << std::flush; };
		std::vector<IxnParamTraj*> bf_ips;
		for (auto d: dims) {
			// Find the basis func dimensions
			bf_ips.clear();
			for (auto bfd: d.basis_func_dims) {
				bf_ips.push_back(_find_ixn_param_by_name(bfd));
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
					bf_ips.push_back(_find_ixn_param_by_name(bfd));
				};
				// Find the basis func
				bf_ptr = _find_basis_func("F_"+denom.name);
				// Find the interaction param
				ixn_param_ptr = _find_ixn_param_by_name(num.name);
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
					bf_ips.push_back(_find_ixn_param_by_name(bfd));
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

	OptProblem::Impl::Impl(Impl&& other) : _time(other._time)
	{
		_copy(other);
		other._reset();
	};

	OptProblem::Impl& OptProblem::Impl::Impl::operator=(Impl&& other)
	{
		if (this != &other)
		{
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;		
	};

	OptProblem::Impl::~Impl() {
		_clean_up();
	};

	/********************
	Helpers for constructors
	********************/

	void OptProblem::Impl::_clean_up()
	{
		// Nothing...
	};
	void OptProblem::Impl::_reset()
	{
		_n_param = 0;
		_ixn_params.clear();
		_bfs.clear();
		_var_terms.clear();
		_hidden_layer_exists = false;
		_hidden_units.clear();
		_species.clear();
		_n_t_soln = 0;
		_t_opt = 0;
		_box_length = 0;
	};

	void OptProblem::Impl::_copy(const Impl& other)
	{
		_n_param = other._n_param;
		_ixn_params = other._ixn_params;
		_bfs = other._bfs;
		_var_terms = other._var_terms;
		_hidden_layer_exists = other._hidden_layer_exists;
		_hidden_units = other._hidden_units;
		_time = other._time;
		_species = other._species;
		_n_t_soln = other._n_t_soln;
		_t_opt = other._t_opt;
		_box_length = other._box_length;
		_latt = other._latt;
	};

	/********************
	Set properties
	********************/

	void OptProblem::Impl::add_hidden_unit(std::vector<std::vector<int>> lattice_idxs, std::string species) 
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
	void OptProblem::Impl::add_hidden_unit(std::vector<int> lattice_idxs, std::string species) 
	{
		// Find sites indicated by connections
		std::vector<Site*> conns;
		for (auto c: lattice_idxs) {
			conns.push_back(_latt.get_site(c));
		};

		// Add
		_add_hidden_unit(conns,species);
	};
	void OptProblem::Impl::_add_hidden_unit(std::vector<Site*> conns, std::string species)
	{
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

		// See if a bias exists for hidden units with this species
		ip = _find_ixn_param_b_by_species(species,false); // can fail
		if (ip) {
			// Tell the bias that this hidden unit exists
			ip->add_hidden_unit(&_hidden_units.back());

			// Tell the hidden unit that this is it's bias
			_hidden_units.back().set_bias(ip);

			// Tell the hidden unit what time it is
			_hidden_units.back().set_t_opt_ptr(&_t_opt);
		};
	};

	/********************
	Validate setup
	********************/

	void OptProblem::Impl::validate_setup() const {
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

	void OptProblem::Impl::solve_ixn_param_traj() {
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

	void OptProblem::Impl::solve_var_traj() {
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

	void OptProblem::Impl::solve(std::vector<std::string> fnames, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolve options)
	{
		// Clear/make directories if needed
		if (options.write) {
			system(("rm -rf " + options.dir_write).c_str());
			system(("mkdir " + options.dir_write).c_str());
			system(("mkdir " + options.dir_write + "/F").c_str());
			system(("mkdir " + options.dir_write + "/ixn_params").c_str());
			system(("mkdir " + options.dir_write + "/moments").c_str());
		};

		// Write the grids
		if (options.write) {
			write_bf_grids(options.dir_write);
			write_t_grid(options.dir_write);
		};

		// For filenames
		std::vector<std::string> fnames_to_use, fnames_remaining;
		std::vector<std::string>::iterator itf;

		// Iterate over optimization steps
		for (int i_opt=0; i_opt<n_opt; i_opt++)
		{
			std::cout << "Opt step " << i_opt << " / " << n_opt-1 << std::endl;

			/*****
			Step 0 - Check nesterov
			*****/

			if (options.nesterov) {
				// If first opt step, do nothing, but set the "prev" point to the current to initialize
				if (i_opt == 0) {
					for (auto itbf=_bfs.begin(); itbf!=_bfs.end(); itbf++) {
						itbf->nesterov_set_prev_equal_curr();
					};
				} else {
					// Move to the intermediate point
					for (auto itbf=_bfs.begin(); itbf!=_bfs.end(); itbf++) {
						itbf->nesterov_move_to_intermediate_pt(i_opt);
					};
				};
			};

			// Write the basis funcs
			if (options.write && !options.write_bf_only_final) {
				write_bfs(options.dir_write+"F/",i_opt);
			};

			/*****
			Step 1 - Solve the current trajectory
			*****/

			if (DIAG_SOLVE) { std::cout << "Solving ixn param" << std::endl; };

			solve_ixn_param_traj();

			// Write
			if (options.write) {
				write_ixn_params(options.dir_write+"ixn_params/",i_opt);
			};

			if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

			/*****
			Step 2 - Solve the variational problem traj
			*****/

			if (DIAG_SOLVE) { std::cout << "Solving var term" << std::endl; };

			solve_var_traj();

			// Write
			/*
			if (options.write) {
				write_var_terms(options.dir_write+"var_terms/",i_opt);
			};
			*/
			
			if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

			/*****
			Step 3 - Pick a random batch
			*****/

			if (DIAG_SOLVE) { std::cout << "Random batch" << std::endl; };

			fnames_to_use.clear();
			fnames_remaining=fnames;
			if (options.use_same_lattice_in_batch) {
				// Use a single lattice - pick a rand
				itf = fnames_remaining.begin();
				std::advance(itf,randI(0,fnames_remaining.size()-1));
				// Add
				fnames_to_use.push_back(*itf);
				for (int i_batch=1; i_batch<batch_size; i_batch++) {
					// Only this one
					fnames_to_use.push_back(fnames_to_use.back());
				};
			} else {
				while (fnames_to_use.size() < batch_size) {
					// Choose
					itf = fnames_remaining.begin();
					std::advance(itf,randI(0,fnames_remaining.size()-1));
					// Add
					fnames_to_use.push_back(*itf);
					// Don't choose again
					fnames_remaining.erase(itf);
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
				if (options.verbose) {
					std::cout << "time: " << _t_opt << std::flush;
				};

				/*****
				Step 5.1 - loop over all samples in the batch
				*****/

				if (DIAG_SOLVE) { std::cout << "   Looping over batch" << std::endl; };

				for (int i_batch=0; i_batch<batch_size; i_batch++) 
				{
					if (options.verbose) {
						std::cout << "." << std::flush;
					};

					/*****
					Step 5.1.1 - Read in sample at this timestep
					*****/

					if (DIAG_SOLVE) { std::cout << "      Read in batch" << std::endl; };

					if (options.awake_visible_are_binary) {
						// Binary
						_latt.read_from_file(fnames_to_use[i_batch] + pad_str(options.time_idx_start+_t_opt,4) + ".txt");
					} else {
						// Probabilistic
						std::cerr << "Error! Probabilistic awake visible units not supported yet!" << std::endl;
						exit(EXIT_FAILURE);
					};

					/*****
					Step 5.1.2 - Activate the hidden units if needed
					*****/

					if (DIAG_SOLVE) { std::cout << "      Activating hidden units" << std::endl; };

					if (_hidden_layer_exists) {
						for (auto ithu = _hidden_units.begin(); ithu != _hidden_units.end(); ithu++) {
							ithu->activate(options.awake_hidden_are_binary);
						};
					};

					/*****
					Step 5.1.3 - Record the awake moments
					*****/

					if (DIAG_SOLVE) { std::cout << "      Record awake moments" << std::endl; };

					for (auto itp = _ixn_params.begin(); itp != _ixn_params.end(); itp++) {
						itp->moments_retrieve_at_time(IxnParamTraj::AWAKE,_t_opt,batch_size);
					};

					/*****
					Step 5.1.4 - CD steps
					*****/

					for (int cd_step=0; cd_step<n_cd_steps; cd_step++)
					{
						// Sample

						if (DIAG_SOLVE) { std::cout << "      Sample" << std::endl; };

						if (cd_step != n_cd_steps-1) {
							if (options.asleep_visible_are_binary) {
								_latt.sample();
							} else {
								std::cerr << "Error! Probabilistic asleep visible units not supported yet!" << std::endl;
								exit(EXIT_FAILURE);
							};
						} else {
							if (options.asleep_final_visible_are_binary) {
								_latt.sample();
							} else {
								std::cerr << "Error! Probabilistic asleep visible units not supported yet!" << std::endl;
								exit(EXIT_FAILURE);
							};
						};

						// Activate the hidden units if needed

						if (DIAG_SOLVE) { std::cout << "      Activating hidden units" << std::endl; };

						if (_hidden_layer_exists) {
							for (auto ithu = _hidden_units.begin(); ithu != _hidden_units.end(); ithu++) {
								if (cd_step != n_cd_steps-1) {
									ithu->activate(options.asleep_hidden_are_binary);
								} else {
									ithu->activate(options.asleep_final_hidden_are_binary);
								};
							};
						};
					};

					/*****
					Step 5.1.5 - Record the asleep moments
					*****/

					if (DIAG_SOLVE) { std::cout << "      Record asleep moments" << std::endl; };

					for (auto itp = _ixn_params.begin(); itp != _ixn_params.end(); itp++) {
						itp->moments_retrieve_at_time(IxnParamTraj::ASLEEP,_t_opt,batch_size);
					};
				};

				if (options.verbose) {
					std::cout << std::endl;
				};

				if (DIAG_SOLVE) { std::cout << "   OK" << std::endl; };
			};

			if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

			// Write the moments
			if (options.write) {
				write_moments(options.dir_write+"moments/",i_opt);
			};

			/*****
			Step 6 - Update the basis funcs
			*****/

			for (auto itbf=_bfs.begin(); itbf!=_bfs.end(); itbf++) {
				itbf->update(_n_t_soln, _time.delta(), dopt);
			};
		};

		// Write the basis funcs one last time
		if (options.write) {
			write_bfs(options.dir_write+"F/",n_opt);
		};
	};



	/********************
	Solve over varying initial conditions
	********************/

	void OptProblem::Impl::solve_varying_ic(std::vector<std::pair<std::string,int>> fnames, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolve options)
	{
		// Clear/make directories if needed
		if (options.write) {
			system(("rm -rf " + options.dir_write).c_str());
			system(("mkdir " + options.dir_write).c_str());
			system(("mkdir " + options.dir_write + "/F").c_str());
			system(("mkdir " + options.dir_write + "/ixn_params").c_str());
			system(("mkdir " + options.dir_write + "/moments").c_str());
		};

		// Write the grids
		if (options.write) {
			write_bf_grids(options.dir_write);
			write_t_grid(options.dir_write);
		};

		// For choosing the batch
		int i_chosen;
		std::vector<std::pair<std::string,int>> fnames_to_use, fnames_remaining;
		std::vector<std::pair<std::string,int>>::iterator itf;

		// Iterate over optimization steps
		for (int i_opt=0; i_opt<n_opt; i_opt++)
		{
			std::cout << "Opt step " << i_opt << " / " << n_opt-1 << std::endl;

			/*****
			Step 0 - Check nesterov
			*****/

			if (options.nesterov) {
				// If first opt step, do nothing, but set the "prev" point to the current to initialize
				if (i_opt == 0) {
					for (auto itbf=_bfs.begin(); itbf!=_bfs.end(); itbf++) {
						itbf->nesterov_set_prev_equal_curr();
					};
				} else {
					// Move to the intermediate point
					for (auto itbf=_bfs.begin(); itbf!=_bfs.end(); itbf++) {
						itbf->nesterov_move_to_intermediate_pt(i_opt);
					};
				};
			};

			// Write the basis funcs
			if (options.write && !options.write_bf_only_final) {
				write_bfs(options.dir_write+"F/",i_opt);
			};

			/*****
			Step 1 - Pick a random batch
			*****/

			if (DIAG_SOLVE) { std::cout << "Random batch" << std::endl; };

			// Clear
			fnames_to_use.clear();
			// Reset
			fnames_remaining=fnames;

			// Files to always use in every batch
			if (options.fnames_used_in_every_batch.size() > 0) {
				fnames_to_use = options.fnames_used_in_every_batch;
			};

			// Go through the batch size
			while (fnames_to_use.size() < batch_size) {
				// Grab a filename
				itf = fnames_remaining.begin();
				i_chosen = randI(0,fnames_remaining.size()-1);
				std::advance(itf,i_chosen);

				// Add that this is to be used
				fnames_to_use.push_back(*itf);

				// Don't choose this again
				fnames_remaining.erase(itf);
			};

			if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

			/*****
			Step 2 - Iterate over batch
			*****/

			if (DIAG_SOLVE) { std::cout << "   Looping over batch" << std::endl; };

			for (int i_batch=0; i_batch<batch_size; i_batch++) 
			{
				if (options.verbose) {
					std::cout << "Doing sample: " << i_batch << " / " << batch_size << " file: " << fnames_to_use[i_batch].first << std::endl;
				};

				/*****
				Step 2.1 - Read the IC
				*****/

				read_init_cond(fnames_to_use[i_batch].first+"../");

				/*****
				Step 2.2 - Solve the current trajectory
				*****/

				if (DIAG_SOLVE) { std::cout << "Solving ixn param" << std::endl; };

				solve_ixn_param_traj();

				// Write
				write_ixn_params(options.dir_write+"ixn_params/",i_opt,fnames_to_use[i_batch].second);

				if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

				/*****
				Step 2.3 - Solve the variational problem traj
				*****/

				if (DIAG_SOLVE) { std::cout << "Solving var term" << std::endl; };

				solve_var_traj();

				// Write
				//write_var_terms(options.dir_write+"var_terms/",i_opt);
				
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
					if (options.verbose) {
						std::cout << "." << std::flush;
					};

					/*****
					Step 2.5.1 - Read in sample at this timestep
					*****/

					if (DIAG_SOLVE) { std::cout << "      Read in batch" << std::endl; };

					if (options.awake_visible_are_binary) {
						// Binary
						_latt.read_from_file(fnames_to_use[i_batch].first + pad_str(options.time_idx_start+_t_opt,4) + ".txt");
					} else {
						// Probabilistic
						std::cerr << "Error! Probabilistic awake visible units are not supported yet." << std::endl;
						exit(EXIT_FAILURE);
					};

					/*****
					Step 2.5.2 - Activate the hidden units if needed
					*****/

					if (DIAG_SOLVE) { std::cout << "      Activating hidden units" << std::endl; };

					if (_hidden_layer_exists) {
						for (auto ithu = _hidden_units.begin(); ithu != _hidden_units.end(); ithu++) {
							ithu->activate(options.awake_hidden_are_binary);
						};
					};

					/*****
					Step 2.5.3 - Record the awake moments
					*****/

					if (DIAG_SOLVE) { std::cout << "      Record awake moments" << std::endl; };

					for (auto itp = _ixn_params.begin(); itp != _ixn_params.end(); itp++) {
						itp->moments_retrieve_at_time(IxnParamTraj::AWAKE,_t_opt);
					};

					/*****
					Step 2.5.4 - CD steps - alternate sampling and activating
					*****/

					for (int cd_step=0; cd_step<n_cd_steps; cd_step++) {

						// Sample

						if (DIAG_SOLVE) { std::cout << "      Sample" << std::endl; };

						if (cd_step != n_cd_steps-1) {
							if (options.asleep_visible_are_binary) {
								_latt.sample();
							} else {
								std::cerr << "Error! Probabilistic asleep visible units not supported yet!" << std::endl;
								exit(EXIT_FAILURE);
							};
						} else {
							if (options.asleep_final_visible_are_binary) {
								_latt.sample();
							} else {
								std::cerr << "Error! Probabilistic asleep visible units not supported yet!" << std::endl;
								exit(EXIT_FAILURE);
							};
						};

						// Activate the hidden units if needed

						if (DIAG_SOLVE) { std::cout << "      Activating hidden units" << std::endl; };

						if (_hidden_layer_exists) {
							for (auto ithu = _hidden_units.begin(); ithu != _hidden_units.end(); ithu++) {
								if (cd_step != n_cd_steps-1) {
									ithu->activate(options.asleep_hidden_are_binary);
								} else {
									ithu->activate(options.asleep_final_hidden_are_binary);
								};
							};
						};

					};

					// Record the asleep moments

					if (DIAG_SOLVE) { std::cout << "      Record asleep moments" << std::endl; };

					for (auto itp = _ixn_params.begin(); itp != _ixn_params.end(); itp++) {
						itp->moments_retrieve_at_time(IxnParamTraj::ASLEEP,_t_opt);
					};

				};

				if (options.verbose) {
					std::cout << std::endl;
				};

				if (DIAG_SOLVE) { std::cout << "   OK" << std::endl; };

				/*****
				Step 2.6 - Write the moments
				*****/

				if (options.write) {
					write_moments(options.dir_write+"moments/",i_opt,fnames_to_use[i_batch].second);
				};

				/*****
				Step 2.7 - Gather the update (but dont commit)
				*****/

				for (auto itbf=_bfs.begin(); itbf!=_bfs.end(); itbf++) {
					itbf->update_gather(_n_t_soln, _time.delta(), dopt);
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
		if (options.write) {
			write_bfs(options.dir_write+"F/",n_opt);
		};
	};


	/********************
	Read some initial conditions
	********************/

	void OptProblem::Impl::read_init_cond(std::string dir) {

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
					ip = _find_ixn_param_by_name(sname);
					ip->set_init_cond(std::stod(sval));

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

	void OptProblem::Impl::write_bf_grids(std::string dir) const {
		for (auto it=_bfs.begin(); it!=_bfs.end(); it++) {
			it->write_grid(dir+"grid_"+it->name()+".txt");
		};
	};
	void OptProblem::Impl::write_t_grid(std::string dir) const {
		_time.write_grid(dir+"grid_time.txt");
	};

	void OptProblem::Impl::write_ixn_params(std::string dir, int idx) const {
		for (auto it = _ixn_params.begin(); it!=_ixn_params.end(); it++) {
			it->write_vals(dir,idx,_n_t_soln);
		};
	};
	void OptProblem::Impl::write_ixn_params(std::string dir, int idx1, int idx2) const {
		for (auto it = _ixn_params.begin(); it!=_ixn_params.end(); it++) {
			it->write_vals(dir,idx1,idx2,_n_t_soln);
		};
	};
	void OptProblem::Impl::write_bfs(std::string dir, int idx) const {
		for (auto it=_bfs.begin(); it!=_bfs.end(); it++) {
			it->write_vals(dir, idx);
		};
	};
	void OptProblem::Impl::write_var_terms(std::string dir, int idx) const {
		for (auto it=_var_terms.begin(); it!=_var_terms.end(); it++) {
			it->write_vals(dir,idx);
		};
	};
	void OptProblem::Impl::write_moments(std::string dir, int idx) const {
		for (auto it = _ixn_params.begin(); it!=_ixn_params.end(); it++) {
			it->write_moments(dir,idx,_n_t_soln);
		};
	};
	void OptProblem::Impl::write_moments(std::string dir, int idx1, int idx2) const {
		for (auto it = _ixn_params.begin(); it!=_ixn_params.end(); it++) {
			it->write_moments(dir,idx1,idx2,_n_t_soln);
		};
	};

	/********************
	Read
	********************/

	void OptProblem::Impl::read_basis_func(std::string bf_name, std::string fname) 
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
	OptProblem - IMPLEMENTATION - PRIVATE
	****************************************/

	/********************
	Search functions
	********************/

	Species* OptProblem::Impl::_find_species(std::string name, bool enforce_success) {
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
	IxnParamTraj* OptProblem::Impl::_find_ixn_param_by_name(std::string name, bool enforce_success) {
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

	IxnParamTraj* OptProblem::Impl::_find_ixn_param_b_by_species(std::string species_name, bool enforce_success) {
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			if (it->is_b_with_species(species_name)) {
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
	IxnParamTraj* OptProblem::Impl::_find_ixn_param_j_by_species(std::string species_name_1, std::string species_name_2, bool enforce_success) {
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
	IxnParamTraj* OptProblem::Impl::_find_ixn_param_w_by_species(std::string species_name, bool enforce_success) {
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
	BasisFunc* OptProblem::Impl::_find_basis_func(std::string name, bool enforce_success) {
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
	VarTermTraj* OptProblem::Impl::_find_var_term(std::string name, bool enforce_success) {
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







	/****************************************
	OptProblem IMPL forwards
	****************************************/

	// Constructor
	OptProblem::OptProblem(std::vector<Dim> dims, double t_max, int n_t, int box_length, int lattice_dim) : _impl(new Impl(dims,t_max,n_t,box_length,lattice_dim)) {};
	OptProblem::OptProblem(OptProblem&& other) = default; // movable but no copies
    OptProblem& OptProblem::operator=(OptProblem&& other) = default; // movable but no copies
	OptProblem::~OptProblem() = default;

	void OptProblem::add_hidden_unit(std::vector<std::vector<int>> lattice_idxs, std::string species) {
		_impl->add_hidden_unit(lattice_idxs,species);
	};
	void OptProblem::add_hidden_unit(std::vector<int> lattice_idxs, std::string species) {
		_impl->add_hidden_unit(lattice_idxs,species);
	};

	void OptProblem::validate_setup() const {
		_impl->validate_setup();
	};

	void OptProblem::solve(std::vector<std::string> fnames, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolve options) {
		_impl->solve(fnames,n_opt,batch_size,n_cd_steps,dopt,options);
	};
	void OptProblem::solve_varying_ic(std::vector<std::pair<std::string,int>> fnames, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolve options) {
		_impl->solve_varying_ic(fnames,n_opt,batch_size,n_cd_steps,dopt,options);
	};

	void OptProblem::read_basis_func(std::string bf_name, std::string fname) {
		_impl->read_basis_func(bf_name,fname);
	};

	void OptProblem::write_bf_grids(std::string dir) const {
		_impl->write_bf_grids(dir);
	};
	void OptProblem::write_t_grid(std::string dir) const {
		_impl->write_t_grid(dir);
	};

	void OptProblem::write_ixn_params(std::string dir, int idx) const {
		_impl->write_ixn_params(dir,idx);
	};
	void OptProblem::write_ixn_params(std::string dir, int idx1, int idx2) const {
		_impl->write_ixn_params(dir,idx1,idx2);
	};
	void OptProblem::write_bfs(std::string dir, int idx) const {
		_impl->write_bfs(dir,idx);
	};
	void OptProblem::write_var_terms(std::string dir, int idx) const {
		_impl->write_var_terms(dir,idx);
	};
	void OptProblem::write_moments(std::string dir, int idx) const {
		_impl->write_moments(dir,idx);
	};
	void OptProblem::write_moments(std::string dir, int idx1, int idx2) const {
		_impl->write_moments(dir,idx1,idx2);
	};

};





