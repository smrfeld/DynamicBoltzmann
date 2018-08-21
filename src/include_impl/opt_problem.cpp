#include "../../include/dynamicboltz_bits/opt_problem.hpp"

// Other headers
#include "../../include/dynamicboltz_bits/dim.hpp"
#include "../../include/dynamicboltz_bits/general.hpp"
#include "../basis_func.hpp"
#include "../counter.hpp"
#include "../hidden_unit.hpp"
#include "../ixn_param_traj.hpp"
#include "../lattice.hpp"
#include "../species.hpp"
#include "../var_term_traj.hpp"

#include <list>
#include <iostream>
#include <ctime>
#include <fstream>
#include <algorithm>
#include <set>
#include <stdlib.h>
#include <sstream>

#define DIAG_SETUP 0
#define DIAG_SOLVE 0
#define DIAG_TIME_SOLVE 0
#define DIAG_SOLVE_RAND_DIV 0

/************************************
* Namespace for dboltz
************************************/

namespace dboltz {

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
		std::list<ConnectionVH> _conn_vh;

		// Counters
		std::list<Counter> _counters;

		// Time dimension
		Grid _time;

		// Species present
		std::list<Species> _species;
		std::list<HiddenSpecies> _hidden_species;

		// Number of steps in this nu solution
		int _n_t_soln;

		// The current time in the optimization step
		int _t_opt;

		// Lattice to hold the current sample of the batch
		Lattice _latt;

		// Add a hidden unit
		void _add_hidden_unit(std::vector<std::string> species_possible, std::vector<Site*> conn_sites, std::vector<std::string> w_params, std::vector<std::string> b_params);
		std::vector<Site*> _get_sites(std::vector<int> &lattice_idxs);
		std::vector<Site*> _get_sites(std::vector<std::vector<int>> &lattice_idxs);

		// Search functions
		Species* _not_nullptr(Species* ptr);
		HiddenSpecies* _not_nullptr(HiddenSpecies *ptr);
		IxnParamTraj* _not_nullptr(IxnParamTraj* ptr);
		BasisFunc* _not_nullptr(BasisFunc* ptr);
		VarTermTraj* _not_nullptr(VarTermTraj* ptr);
		Counter* _not_nullptr(Counter* ptr);
		// Find species
		Species* _find_species(std::string name);
		// Find hidden species
		HiddenSpecies* _find_hidden_species(std::string name);
		// Find ixn param
		IxnParamTraj* _find_ixn_param(std::string name);
		// Find basis function
		BasisFunc* _find_basis_func(std::string name);
		// Find var term traj
		VarTermTraj* _find_var_term(std::string name);
		// Find counter
		Counter* _find_ctr_by_species(std::string s);
		Counter* _find_ctr_by_species(std::string s1, std::string s2);
		Counter* _find_ctr_by_species(std::string s1, std::string s2, std::string s3);
		Counter* _find_ctr_by_species(std::string s1, std::string s2, std::string s3, std::string s4);

		// Add a counter for some species or nns
		void _add_counter(std::string s);
		void _add_counter(std::string s1, std::string s2);
		void _add_counter(std::string s1, std::string s2, std::string s3);
		void _add_counter(std::string s1, std::string s2, std::string s3, std::string s4);

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const Impl& other);

	public:

		/********************
		Constructor
		********************/

		Impl(std::vector<Dim> dims, std::vector<std::string> species_visible, std::vector<std::string> species_hidden, double t_max, int n_t, int box_length, int lattice_dim);
		Impl(const Impl& other);
		Impl(Impl&& other);
	    Impl& operator=(const Impl& other);
	    Impl& operator=(Impl&& other);
		~Impl();

		/********************
		Change the time limit
		********************/

		void set_n_t(int n_t);

		/********************
		Set IC for ixn param
		********************/

		void set_ic_for_ixn_param(std::string param_name, double val);

		/********************
		Set a fixed awake moment
		********************/

		void set_fixed_awake_moment_for_dim(std::string param_name, std::vector<double> vals);

		/********************
		Set properties	
		********************/

		// Any dim
		void add_hidden_unit(std::vector<std::string> species_possible, std::vector<std::vector<int>> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params);
		// 1D specific
		void add_hidden_unit(std::vector<std::string> species_possible, std::vector<int> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params);

		/********************
		Validate setup by printing
		********************/

		void validate_setup() const;
		void validate_graph() const;

		/********************
		Solve interaction parameter traj
		********************/

		void solve_ixn_param_traj(int t_end=0);

		/********************
		Solve for variational trajectory
		********************/

		void solve_var_traj();
		void solve_var_traj_from_zero(int it_start, int it_end);

		/********************
		Solve
		********************/

		void solve(std::vector<std::string> fnames, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolve options);
		void solve(std::vector<std::vector<std::string>> fnames, int n_opt, int n_cd_steps, double dopt, OptionsSolve options);
		void _solve(std::vector<std::vector<std::string>> fnames_to_use, int n_opt, int n_cd_steps, double dopt, OptionsSolve options);

		void solve_varying_ic(std::vector<FName> fnames, std::vector<FName> fnames_used_in_every_batch, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolve options);
		void solve_varying_ic(std::vector<std::vector<FName>> fnames, int n_opt, int n_cd_steps, double dopt, OptionsSolve options);
		void _solve_varying_ic(std::vector<std::vector<FName>> fnames_to_use, int n_opt, int n_cd_steps, double dopt, OptionsSolve options);

		void solve_rand_div(std::vector<std::string> fname_collection, int n_opt, int batch_size, int n_cd_steps, double dopt, int n_divisions, OptionsSolve options);
		void solve_rand_div(std::vector<std::vector<std::string>> fnames, int n_opt, int n_cd_steps, double dopt, int n_divisions, OptionsSolve options);
		void _solve_rand_div(std::vector<std::vector<std::string>> fnames_to_use, int n_opt, int n_cd_steps, double dopt, int n_divisions, OptionsSolve options);

		/********************
		Read basis func
		********************/

		void read_basis_func(std::string bf_name, std::string fname);

		void read_init_cond_for_ixn_param(std::string ixn_func_name, std::string fname, int line_idx=-1);

		/********************
		Read some initial conditions
		********************/

		void read_init_cond(std::string dir);

		/********************
		Write
		********************/

		void write_bf_grids(std::string dir) const;
		void write_t_grid(std::string dir) const;

		void write_ixn_params(std::string dir, int idx_opt_step) const;
		void write_ixn_params(std::string dir, int idx_opt_step, std::vector<int> idxs) const;
		void write_bfs(std::string dir, int idx_opt_step) const;
		void write_var_terms(std::string dir, int idx_opt_step) const;
		void write_moments(std::string dir, int idx_opt_step) const;
		void write_moments(std::string dir, int idx_opt_step, std::vector<int> idxs) const;
	};








































	/****************************************
	OptProblem - IMPLEMENTATION DEFINITIONS
	****************************************/

	/********************
	Constructor
	********************/

	OptProblem::Impl::Impl(std::vector<Dim> dims, std::vector<std::string> species_visible, std::vector<std::string> species_hidden, double t_max, int n_t, int box_length, int lattice_dim) : _latt(lattice_dim,box_length), _time("time",0.0,t_max,n_t)
	{
		// Set parameters
		if (DIAG_SETUP) { std::cout << "Copying params..." << std::flush; };
		_n_param = dims.size();
		_t_opt = 0;
		_n_t_soln = 0;
		_hidden_layer_exists = false;
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Fix dims that have any species specification
		if (DIAG_SETUP) { std::cout << "Fixing incomplete dim specifications (for any species)..." << std::flush; };
		for (auto &d: dims) {
			if (d.any_species()) {
				if (d.type()==H) {
					for (auto s: species_visible) {
						d.add_species_h(s);
					};
				} else if (d.type()==B) {
					for (auto s: species_hidden) {
						d.add_species_b(s);
					};
				} else if (d.type()==J) {
					for (auto s1: species_visible) {
						for (auto s2: species_visible) {
							d.add_species_J(s1,s2);
						};
					};
				} else if (d.type()==W) {
					for (auto sv: species_visible) {
						for (auto sh: species_hidden) {
							d.add_species_W(sv,sh);
						};
					};
				};
			};
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Check if a hidden layer exists
		if (DIAG_SETUP) { std::cout << "Checking if hidden layer exists..." << std::flush; };
		for (auto const &d: dims) {
			if (d.type()==W || d.type()==B) {
				_hidden_layer_exists = true;
				break;
			};
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Create the visible and hidden species
		if (DIAG_SETUP) { std::cout << "Create species..." << std::flush; };
		for (auto const &s: species_visible) {
			_species.push_back(Species(s));

			// Add to the lattice
			_latt.add_species_possibility(&(_species.back()));
		};
		for (auto const &s: species_hidden) {
			_hidden_species.push_back(HiddenSpecies(s));
		};		
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Create the interaction params
		if (DIAG_SETUP) { std::cout << "Create ixn params..." << std::flush; };
		Species *sp,*sp1,*sp2,*sp3;
		HiddenSpecies *sph;
		for (auto const &d: dims) {

			if (d.type()==H) {

				// Create
				_ixn_params.push_back(IxnParamTraj(d.name(),Hp,d.min(),d.max(),d.n(),d.init(),n_t,&_t_opt));

				for (auto s: d.get_species_h()) {
					// Find and add the species
					sp = _not_nullptr(_find_species(s));
					_ixn_params.back().add_species(sp);
					// Add ixn to the species
					sp->add_h_ptr(&_ixn_params.back());
					// Make sure we have a counter
					_add_counter(s);					
				};

			} else if (d.type()==J) { 

				// Create
				_ixn_params.push_back(IxnParamTraj(d.name(),Jp,d.min(),d.max(),d.n(),d.init(),n_t,&_t_opt));

				for (auto s: d.get_species_J()) {
					// Find and add the species
					sp1 = _not_nullptr(_find_species(s[0]));
					sp2 = _not_nullptr(_find_species(s[1]));
					_ixn_params.back().add_species(sp1,sp2);
					// Add ixn to the species
					sp1->add_j_ptr(sp2,&_ixn_params.back());
					sp2->add_j_ptr(sp1,&_ixn_params.back());
					// Make sure we have a counter
					_add_counter(s[0],s[1]);					
				};

				// Init structure for lattice
				_latt.init_nn_structure();

			} else if (d.type()==K) {

				// Check that Lattice dim is 1 - only 1 is allowed!
				if (lattice_dim != 1) {
					std::cerr << "Error: triplets are currently only supported for lattice of dim 1" << std::endl;
					exit(EXIT_FAILURE);
				};

				// Create
				_ixn_params.push_back(IxnParamTraj(d.name(),Kp,d.min(),d.max(),d.n(),d.init(),n_t,&_t_opt));

				for (auto s: d.get_species_K()) {
					// Find and add the species
					sp1 = _not_nullptr(_find_species(s[0]));
					sp2 = _not_nullptr(_find_species(s[1]));
					sp3 = _not_nullptr(_find_species(s[2]));
					_ixn_params.back().add_species(sp1,sp2,sp3);
					// Add ixn to the species
					sp1->add_k_ptr(sp2,sp3,&_ixn_params.back());
					sp2->add_k_ptr(sp1,sp3,&_ixn_params.back());
					sp3->add_k_ptr(sp1,sp2,&_ixn_params.back());
					// Make sure we have a counter
					_add_counter(s[0],s[1],s[2]);			
				};

				// Init structure for lattice
				_latt.init_triplet_structure();

			} else if (d.type()==W) { 

				// Create
				_ixn_params.push_back(IxnParamTraj(d.name(),Wp,d.min(),d.max(),d.n(),d.init(),n_t,&_t_opt));

				for (auto s: d.get_species_W()) {
					// Find and add the species
					sp = _not_nullptr(_find_species(s[0]));
					sph = _not_nullptr(_find_hidden_species(s[1]));
					_ixn_params.back().add_species(sp,sph);

					// No need to add to species - handled by ConnectionVH class
				};

			} else if (d.type()==B) {

				// Create
				_ixn_params.push_back(IxnParamTraj(d.name(),Bp,d.min(),d.max(),d.n(),d.init(),n_t,&_t_opt));

				for (auto s: d.get_species_b()) {
					// Find and add the species
					sph = _not_nullptr(_find_hidden_species(s));
					_ixn_params.back().add_species(sph);

					// No need to add to species - handled by ConnectionVH class
				};

			};
		};

		// Create the basis functions
		if (DIAG_SETUP) { std::cout << "Create basis funcs..." << std::flush; };
		std::vector<IxnParamTraj*> bf_ips;
		IxnParamTraj *iptr;
		for (auto const &d: dims) {
			// Find the basis func dimensions
			bf_ips.clear();
			for (auto bfd: d.basis_func_dims()) {
				bf_ips.push_back(_not_nullptr(_find_ixn_param(bfd)));
			};
			// Find the ixn param
			iptr = _not_nullptr(_find_ixn_param(d.name()));
			// Make the basis function
			_bfs.push_back(BasisFunc("F_"+d.name(),iptr,bf_ips));
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Add the basis functions to the interaction params
		if (DIAG_SETUP) { std::cout << "Add basis func to ixn params..." << std::flush; };
		BasisFunc* bf_ptr=nullptr;
		for (auto itp=_ixn_params.begin(); itp!=_ixn_params.end(); itp++) {
			// Find the basis func
			bf_ptr = _not_nullptr(_find_basis_func("F_"+itp->name()));
			// Set
			itp->set_basis_func_ptr(bf_ptr);
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Create the variational terms
		if (DIAG_SETUP) { std::cout << "Create var terms..." << std::flush; };
		BasisFunc *num_bf_ptr = nullptr;
		IxnParamTraj *ixn_param_ptr = nullptr;
		// Go through numerators
		for (auto const &num: dims) {
			// Go through denominators
			for (auto const &denom: dims) {
				// Find the basis func dimensions
				bf_ips.clear();
				for (auto bfd: denom.basis_func_dims()) {
					bf_ips.push_back(_not_nullptr(_find_ixn_param(bfd)));
				};
				// Find the basis func
				bf_ptr = _not_nullptr(_find_basis_func("F_"+denom.name()));
				// Find the interaction param
				ixn_param_ptr = _not_nullptr(_find_ixn_param(num.name()));
				// Find the basis func corresponding to the numerator
				num_bf_ptr = _not_nullptr(_find_basis_func("F_"+num.name()));
				// Create the var term
				_var_terms.push_back(VarTermTraj("var_"+num.name()+"_wrt_F_"+denom.name(), ixn_param_ptr, bf_ptr, bf_ips, num_bf_ptr, num.basis_func_dims().size(), n_t));
			};
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Add the pointers to variational terms needed to update this var term
		if (DIAG_SETUP) { std::cout << "Add ptrs to var term..." << std::flush; };
		VarTermTraj* vt_ptr=nullptr;
		// Go through numerators
		for (auto const &num: dims) {
			// Go through denoms
			for (auto denom=_bfs.begin(); denom!=_bfs.end(); denom++) {
				// Find the ixn params that are arguments to the num's basis func
				bf_ips.clear();
				for (auto bfd: num.basis_func_dims()) {
					bf_ips.push_back(_not_nullptr(_find_ixn_param(bfd)));
				};
				// Find the variational term
				vt_ptr = _not_nullptr(_find_var_term("var_"+num.name()+"_wrt_"+denom->name()));
				// Find the variational terms needed to update this one
				for (auto ip_ptr: bf_ips) {
					vt_ptr->add_update_ptr(_not_nullptr(_find_var_term("var_"+ip_ptr->name()+"_wrt_"+denom->name())));
				};
			};
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Add pointers to the basis functions needed to update them
		if (DIAG_SETUP) { std::cout << "Add ptrs to basis func..." << std::flush; };
		for (auto itbf=_bfs.begin(); itbf!=_bfs.end(); itbf++) {
			for (auto itp=_ixn_params.begin(); itp!=_ixn_params.end(); itp++) {
				// Find the variational term
				vt_ptr = _not_nullptr(_find_var_term("var_"+itp->name()+"_wrt_"+itbf->name()));
				// Add
				itbf->add_update_ptrs(&*itp,vt_ptr);
			};
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

	};
	OptProblem::Impl::Impl(const Impl& other) : _time(other._time), _latt(other._latt)
	{
		_copy(other);
	};
	OptProblem::Impl::Impl(Impl&& other) : _time(other._time), _latt(other._latt)
	{
		_copy(other);
		other._reset();
	};
    OptProblem::Impl& OptProblem::Impl::Impl::operator=(const Impl& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
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
		_conn_vh.clear();
		_species.clear();
		_hidden_species.clear();
		_n_t_soln = 0;
		_t_opt = 0;
	};

	void OptProblem::Impl::_copy(const Impl& other)
	{
		_n_param = other._n_param;
		_ixn_params = other._ixn_params;
		_bfs = other._bfs;
		_var_terms = other._var_terms;
		_hidden_layer_exists = other._hidden_layer_exists;
		_hidden_units = other._hidden_units;
		_conn_vh = other._conn_vh;
		_time = other._time;
		_species = other._species;
		_hidden_species = other._hidden_species;
		_n_t_soln = other._n_t_soln;
		_t_opt = other._t_opt;
		_latt = other._latt;
	};

	/********************
	Change the time limit
	********************/

	void OptProblem::Impl::set_n_t(int n_t) {
		// Change time object
		_time.set_n(n_t);
		// Ixn terms
		for (auto itp=_ixn_params.begin(); itp!=_ixn_params.end(); itp++) {
			itp->set_n_t(n_t);
		};
		// Var terms
		for (auto itv=_var_terms.begin(); itv!=_var_terms.end(); itv++) {
			itv->set_n_t(n_t);
		};
	};

	/********************
	Set IC for ixn param
	********************/

	void OptProblem::Impl::set_ic_for_ixn_param(std::string param_name, double val) {
		// Find
		IxnParamTraj *ip = _not_nullptr(_find_ixn_param(param_name));
		// Set
		ip->set_init_cond(val);
	};

	/********************
	Set a fixed awake moment
	********************/

	void OptProblem::Impl::set_fixed_awake_moment_for_dim(std::string param_name, std::vector<double> vals) {
		// Find
		IxnParamTraj *ip = _not_nullptr(_find_ixn_param(param_name));
		// Set
		ip->set_fixed_awake_moment(vals);	
	};

	/********************
	Find hidden unit connections
	********************/

	// Any dim
	std::vector<Site*> OptProblem::Impl::_get_sites(std::vector<std::vector<int>> &lattice_idxs) {
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
		return conns;
	};

	// 1D specific
	std::vector<Site*> OptProblem::Impl::_get_sites(std::vector<int> &lattice_idxs) {
		std::vector<Site*> conns;
		for (auto c: lattice_idxs) {
			conns.push_back(_latt.get_site(c));
		};
		return conns;
	};

	/********************
	Add hidden unit
	********************/

	void OptProblem::Impl::add_hidden_unit(std::vector<std::string> species_possible, std::vector<std::vector<int>> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params) {
		// Find sites indicated by connections
		std::vector<Site*> conns = _get_sites(lattice_idxs);
		// Make
		_add_hidden_unit(species_possible, conns, w_params, b_params);
	};
	void OptProblem::Impl::add_hidden_unit(std::vector<std::string> species_possible, std::vector<int> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params) {
		// Find sites indicated by connections
		std::vector<Site*> conns = _get_sites(lattice_idxs);
		// Make
		_add_hidden_unit(species_possible, conns, w_params, b_params);
	};

	/********************
	Add hidden unit internal
	********************/

	void OptProblem::Impl::_add_hidden_unit(std::vector<std::string> species_possible, std::vector<Site*> conn_sites, std::vector<std::string> w_params, std::vector<std::string> b_params) {

		// Find the ixn params listed
		std::vector<IxnParamTraj*> ip_w;
		std::vector<IxnParamTraj*> ip_b;
		IxnParamTraj *ip;
		for (auto w: w_params) {
			ip = _not_nullptr(_find_ixn_param(w));
			ip_w.push_back(ip);
		};
		for (auto b: b_params) {
			ip = _not_nullptr(_find_ixn_param(b));
			ip_b.push_back(ip);
		};

		// Make hidden unit
		_hidden_units.push_back(HiddenUnit());
		HiddenUnit* hup = &_hidden_units.back();

		// Add all allowed hidden species to this hidden unit
		HiddenSpecies *hsp;
		for (auto hsp_name: species_possible) {
			hsp = _not_nullptr(_find_hidden_species(hsp_name));
			hup->add_hidden_species_possibility(hsp);
		};

		// Add the biases to the hidden units and vice versa
		for (auto b: ip_b) {
			// bias to hidden unit
			hup->add_bias(b);
			// hidden unit to bias
			b->add_hidden_unit(hup);
		};

		// Make the connections
		ConnectionVH *cvh;
		for (auto s: conn_sites) {
			// Make the connection
			_conn_vh.push_back(ConnectionVH(s,hup,ip_w));
			cvh = &_conn_vh.back();

			// Add the connection to the hidden unit
			hup->add_visible_hidden_conn(cvh);

			// Add the connection to the lattice site
			s->add_visible_hidden_conn(cvh);

			// Add to ixn params
			for (auto w: ip_w) {
				w->add_visible_hidden_connection(s,hup);
			};
		};
	};

	/********************
	Validate setup
	********************/

	void OptProblem::Impl::validate_setup() const {
		std::cout << "------------------------" << std::endl;
		for (auto it: _species) {
			//it.validate_setup();
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
	void OptProblem::Impl::validate_graph() const {
		_latt.validate_graph();
	};

	/********************
	Solve interaction parameter traj
	********************/

	void OptProblem::Impl::solve_ixn_param_traj(int t_end) {
		// Number of time points before the solution has gone out of bounds
		_n_t_soln = 1;

		// Have we gone out of the domain?
		bool in_domain;

		// Max time
		int t_end_use = t_end;
		if (t_end==0) {
			t_end_use = _time.n();	
		};

		// Go through all times
		for (int it=1; it<t_end_use; it++) {
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
				itv->calculate_at_time(it,_time.delta());
			};
		};
	};

	void OptProblem::Impl::solve_var_traj_from_zero(int it_start, int it_end) {
		// Set all to zero at it_start
		for (auto itv=_var_terms.begin(); itv!=_var_terms.end(); itv++) {
			itv->set_to_zero_at_time(it_start);
		};
		// Go through times and solve
		for (int it=it_start+1; it<it_end; it++) {
			// Go through all var terms
			for (auto itv=_var_terms.begin(); itv!=_var_terms.end(); itv++) {
				itv->calculate_at_time(it,_time.delta());
			};
		};
	};










	/****************************************
	Solve -- Main opt loop
	****************************************/

	void OptProblem::Impl::solve(std::vector<std::string> fnames, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolve options) {

		// Pick filenames
		std::vector<std::vector<std::string>> fnames_to_use;
		std::vector<std::string> fnames_batch;
		std::vector<std::string> fnames_remaining;
		std::vector<std::string>::iterator itf;

		for (int i_opt=0; i_opt < n_opt; i_opt++) {

			// Clear batch
			fnames_batch.clear();

			// Reset remaining
			fnames_remaining=fnames;

			// Choose batch
			if (options.use_same_lattice_in_batch) {
				// Use a single lattice - pick a rand
				itf = fnames_remaining.begin();
				std::advance(itf,randI(0,fnames_remaining.size()-1));
				// Add
				fnames_batch.push_back(*itf);
				// Only this one
				for (int i_batch=1; i_batch<batch_size; i_batch++) {
					fnames_batch.push_back(fnames_batch.back());
				};
			} else {
				while (fnames_batch.size() < batch_size) {
					// Choose a random one
					itf = fnames_remaining.begin();
					std::advance(itf,randI(0,fnames_remaining.size()-1));
					// Add
					fnames_batch.push_back(*itf);
					// Don't choose again
					fnames_remaining.erase(itf);
				};
			};

			/*
			for (auto s: fnames_batch) {
				std::cout << s << " ";
			};
			std::cout << "" << std::endl;
			*/
			
			// Add the batch
			fnames_to_use.push_back(fnames_batch);
		};

		// Main solve function
		_solve(fnames_to_use, n_opt, n_cd_steps, dopt, options);
	};

	void OptProblem::Impl::solve(std::vector<std::vector<std::string>> fname_collection, int n_opt, int n_cd_steps, double dopt, OptionsSolve options) {

		// Pick filenames
		std::vector<std::vector<std::string>> fnames_to_use;
		std::vector<std::string> fnames_batch;
		std::vector<std::string>::iterator itf;

		for (int i_opt=0; i_opt < n_opt; i_opt++) {
			// Pick a random one
			fnames_batch = fname_collection[randI(0,fname_collection.size()-1)];

			// Add
			fnames_to_use.push_back(fnames_batch);
		};

		// Main solve function
		_solve(fnames_to_use, n_opt, n_cd_steps, dopt, options);

	};

	void OptProblem::Impl::_solve(std::vector<std::vector<std::string>> fnames_to_use, int n_opt, int n_cd_steps, double dopt, OptionsSolve options)
	{
		// Clear/make directories if needed
		if (options.clear_dir) {
			system(("rm -rf " + options.dir_write).c_str());
		};
		if (options.write) {
			system(("mkdir -p " + options.dir_write).c_str());
			system(("mkdir -p " + options.dir_write + "F").c_str());
			system(("mkdir -p " + options.dir_write + "ixn_params").c_str());
			system(("mkdir -p " + options.dir_write + "moments").c_str());
		};
		if (options.write_var_terms) {
			system(("mkdir -p " + options.dir_write + "var_terms").c_str());
		};

		// Write the grids
		if (options.write) {
			write_bf_grids(options.dir_write);
			write_t_grid(options.dir_write);
		};

		// Opt step with offset
		int i_opt;

		// Check exp decay / cutoff times
		if (options.exp_decay) {
			// Check
			if (options.exp_decay_t0_values.size() != n_opt || options.exp_decay_lambda_values.size() != n_opt) {
				std::cerr << "Error! In exp decay mode, must provide t0 and lambda values for all optimization timesteps." << std::endl;
				exit(EXIT_FAILURE);
			};
		};
		if (options.time_cutoff) {
			// Check
			if (options.time_cutoff_start_values.size() != n_opt || options.time_cutoff_end_values.size() != n_opt) {
				std::cerr << "Error! In time cutoff mode, must provide values for all optimization timesteps." << std::endl;
				exit(EXIT_FAILURE);
			};
		};

		// Convert l2 reg for params
		std::map<IxnParamTraj*,double> l2_lambda_params, l2_reg_centers;
		if (options.l2_reg_params_mode) {
			for (auto const &pr: options.l2_lambda_params) {
				l2_lambda_params[_not_nullptr(_find_ixn_param(pr.first))] = pr.second;
			};
			for (auto const &pr: options.l2_reg_centers) {
				l2_reg_centers[_not_nullptr(_find_ixn_param(pr.first))] = pr.second;
			};
		};

		// Times possibly
		clock_t t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;

		// Values for exp decay
		double exp_decay_t0=0., exp_decay_lambda=0.;

		// Start/stop time for integration
		int t_start,t_end;

		// Iterate over optimization steps
		for (int i_opt_from_zero=0; i_opt_from_zero<n_opt; i_opt_from_zero++)
		{
			std::cout << "Opt step " << i_opt_from_zero << " / " << n_opt-1 << std::endl;

			if (DIAG_TIME_SOLVE) { t0 = clock(); };

			// Offset
			i_opt = i_opt_from_zero + options.opt_idx_start_writing;

			// Exp decay
			if (options.exp_decay) {
				exp_decay_t0 = options.exp_decay_t0_values[i_opt_from_zero];
				exp_decay_lambda = options.exp_decay_lambda_values[i_opt_from_zero];
			};

			// Start and stop time
			if (options.time_cutoff) {
				t_start = options.time_cutoff_start_values[i_opt_from_zero];
				t_end = options.time_cutoff_end_values[i_opt_from_zero];
			} else {
				t_start = 0;
				t_end = _time.n();
			};

			/*****
			Step 0 - Check nesterov
			*****/

			if (DIAG_SOLVE) { std::cout << "Taking Nesterov step" << std::endl; };

			if (options.nesterov) {
				// If first opt step, do nothing, but set the "prev" point to the current to initialize
				if (i_opt_from_zero == 0) {
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

			if (DIAG_TIME_SOLVE) {
				t1 = clock();
				std::cout << "	nesterov: " << double(t1 - t0) / CLOCKS_PER_SEC << std::endl;
			};

			/*****
			Step 1 - Solve the current trajectory
			*****/

			if (DIAG_SOLVE) { std::cout << "Solving ixn param" << std::endl; };

			solve_ixn_param_traj(t_end);

			// Write
			if (options.write && options.write_ixn_params) {
				write_ixn_params(options.dir_write+"ixn_params/",i_opt);
			};

			if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

			if (DIAG_TIME_SOLVE) {
				t2 = clock();
				std::cout << "	traj: " << double(t2 - t1) / CLOCKS_PER_SEC << std::endl;
			};

			/*****
			Step 2 - Solve the variational problem traj
			*****/

			if (DIAG_SOLVE) { std::cout << "Solving var term" << std::endl; };

			solve_var_traj();

			// Write
			if (options.write && options.write_var_terms) {
				write_var_terms(options.dir_write+"var_terms/",i_opt);
			};
			
			if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

			if (DIAG_TIME_SOLVE) {
				t3 = clock();
				std::cout << "	var: " << double(t3 - t2) / CLOCKS_PER_SEC << std::endl;
			};

			/*****
			Step 4 - reset the moments at all times
			*****/

			if (DIAG_SOLVE) { std::cout << "Reset moments" << std::endl; };

			for (auto itp = _ixn_params.begin(); itp != _ixn_params.end(); itp++) {
				itp->moments_reset();
			};

			if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

			if (DIAG_TIME_SOLVE) {
				t4 = clock();
				std::cout << "	reset: " << double(t4 - t3) / CLOCKS_PER_SEC << std::endl;
			};

			/*****
			Step 5 - Go through all times
			*****/

			if (DIAG_SOLVE) { std::cout << "Looping over times" << std::endl; };

			// Use the class variable _t_opt to iterate
			for (_t_opt=t_start; _t_opt < _n_t_soln; _t_opt++)
			{
				if (options.verbose) {
					std::cout << "time: " << _t_opt << std::flush;
				};

				/*****
				Step 5.1 - loop over all samples in the batch
				*****/

				if (DIAG_SOLVE) { std::cout << "   Looping over batch" << std::endl; };

				for (int i_batch=0; i_batch<fnames_to_use[i_opt_from_zero].size(); i_batch++) 
				{
					if (options.verbose) {
						std::cout << "." << std::flush;
					};


					if (DIAG_TIME_SOLVE) { t5 = clock(); };

					/*****
					Step 5.1.1 - Read in sample at this timestep
					*****/

					if (DIAG_SOLVE) { std::cout << "      Read in batch" << std::endl; };

					if (options.awake_visible_are_binary) {
						// Binary
						_latt.read_from_file(fnames_to_use[i_opt_from_zero][i_batch] + pad_str(options.time_idx_start_reading+_t_opt,4) + ".txt");
					} else {
						// Probabilistic
						std::cerr << "Error! Probabilistic awake visible units not supported yet!" << std::endl;
						exit(EXIT_FAILURE);
					};
 	
 					if (DIAG_TIME_SOLVE) {
						t6 = clock();
						std::cout << "	read: " << double(t6 - t5) / CLOCKS_PER_SEC << std::endl;
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

					if (DIAG_TIME_SOLVE) {
						t7 = clock();
						std::cout << "	activated: " << double(t7 - t6) / CLOCKS_PER_SEC << std::endl;
					};

					/*****
					Step 5.1.3 - Record the awake moments
					*****/

					if (DIAG_SOLVE) { std::cout << "      Record awake moments" << std::endl; };

					for (auto itp = _ixn_params.begin(); itp != _ixn_params.end(); itp++) {
						itp->moments_retrieve_at_time(IxnParamTraj::AWAKE,_t_opt,fnames_to_use[i_opt_from_zero].size());
					};

					if (DIAG_TIME_SOLVE) {
						t8 = clock();
						std::cout << "	recorded: " << double(t8 - t7) / CLOCKS_PER_SEC << std::endl;
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
								_latt.sample(_t_opt);
							} else {
								std::cerr << "Error! Probabilistic asleep visible units not supported yet!" << std::endl;
								exit(EXIT_FAILURE);
							};
						} else {
							if (options.asleep_final_visible_are_binary) {
								_latt.sample(_t_opt);
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

					if (DIAG_TIME_SOLVE) {
						t9 = clock();
						std::cout << "	CD: " << double(t9 - t8) / CLOCKS_PER_SEC << std::endl;
					};

					/*****
					Step 5.1.5 - Record the asleep moments
					*****/

					if (DIAG_SOLVE) { std::cout << "      Record asleep moments" << std::endl; };

					for (auto itp = _ixn_params.begin(); itp != _ixn_params.end(); itp++) {
						itp->moments_retrieve_at_time(IxnParamTraj::ASLEEP,_t_opt,fnames_to_use[i_opt_from_zero].size());
					};

					if (DIAG_TIME_SOLVE) {
						t10 = clock();
						std::cout << "	recorded: " << double(t10 - t9) / CLOCKS_PER_SEC << std::endl;
					};
				};

				if (options.verbose) {
					std::cout << std::endl;
				};

				if (DIAG_SOLVE) { std::cout << "   OK" << std::endl; };
			};

			if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

			// Write the moments
			if (options.write && options.write_moments) {
				write_moments(options.dir_write+"moments/",i_opt);
			};

			/*****
			Step 6 - Update the basis funcs
			*****/

			for (auto itbf=_bfs.begin(); itbf!=_bfs.end(); itbf++) {
				itbf->update(t_start, _n_t_soln, _time.delta(), dopt, options.exp_decay, exp_decay_t0, exp_decay_lambda, options.l2_reg_params_mode, l2_lambda_params, l2_reg_centers);
			};
		};

		// Write the basis funcs one last time
		if (options.write) {
			write_bfs(options.dir_write+"F/",n_opt+options.opt_idx_start_writing);
		};
	};












	/****************************************
	Solve over varying initial conditions
	****************************************/

	void OptProblem::Impl::solve_varying_ic(std::vector<FName> fnames, std::vector<FName> fnames_used_in_every_batch, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolve options) {

		// Pick filenames
		std::vector<std::vector<FName>> fnames_to_use;
		std::vector<FName> fnames_batch;
		std::vector<FName> fnames_remaining;
		std::vector<FName>::iterator itf;
		int i_chosen;

		for (int i_opt=0; i_opt<n_opt; i_opt++) {

			// Clear batch
			fnames_batch.clear();

			// Reset
			fnames_remaining=fnames;

			// Files to always use in every batch
			if (fnames_used_in_every_batch.size() > 0) {
				fnames_batch = fnames_used_in_every_batch;
			};

			// Go through the batch size
			while (fnames_batch.size() < batch_size) {
				// Grab a filename
				itf = fnames_remaining.begin();
				i_chosen = randI(0,fnames_remaining.size()-1);
				std::advance(itf,i_chosen);

				// Add that this is to be used
				fnames_batch.push_back(*itf);

				// Don't choose this again
				fnames_remaining.erase(itf);
			};

			// Add the batch
			fnames_to_use.push_back(fnames_batch);
		};

		// Main solve function
		_solve_varying_ic(fnames_to_use, n_opt, n_cd_steps, dopt, options);
	};

	void OptProblem::Impl::solve_varying_ic(std::vector<std::vector<FName>> fname_collection, int n_opt, int n_cd_steps, double dopt, OptionsSolve options) {

		// Pick filenames
		std::vector<std::vector<FName>> fnames_to_use;
		std::vector<FName> fnames_batch;

		for (int i_opt=0; i_opt < n_opt; i_opt++) {
			// Pick a random one
			fnames_batch = fname_collection[randI(0,fname_collection.size()-1)];

			// Add
			fnames_to_use.push_back(fnames_batch);
		};

		// Main solve function
		_solve_varying_ic(fnames_to_use, n_opt, n_cd_steps, dopt, options);
	};

	void OptProblem::Impl::_solve_varying_ic(std::vector<std::vector<FName>> fnames_to_use, int n_opt, int n_cd_steps, double dopt, OptionsSolve options)
	{
		// Clear/make directories if needed
		if (options.clear_dir) {
			system(("rm -rf " + options.dir_write).c_str());
		};
		if (options.write) {
			system(("mkdir -p " + options.dir_write).c_str());
			system(("mkdir -p " + options.dir_write + "F").c_str());
			if (options.write_ixn_params) {
				system(("mkdir -p " + options.dir_write + "ixn_params").c_str());
			};
			if (options.write_moments) {
				system(("mkdir -p " + options.dir_write + "moments").c_str());
			};
			if (options.write_var_terms) {
				system(("mkdir -p " + options.dir_write + "var_terms").c_str());
			};
		};

		// Write the grids
		if (options.write) {
			write_bf_grids(options.dir_write);
			write_t_grid(options.dir_write);
		};

		// Convert l2 reg for params
		std::map<IxnParamTraj*,double> l2_lambda_params, l2_reg_centers;
		if (options.l2_reg_params_mode) {
			for (auto const &pr: options.l2_lambda_params) {
				l2_lambda_params[_not_nullptr(_find_ixn_param(pr.first))] = pr.second;
			};
			for (auto const &pr: options.l2_reg_centers) {
				l2_reg_centers[_not_nullptr(_find_ixn_param(pr.first))] = pr.second;
			};
		};

		// Optimization step translated by the offset
		int i_opt_translated;

		// Exp decay terms
		double exp_decay_t0=0., exp_decay_lambda=0.;

		// Start/stop time
		int t_start, t_end;

		// Split idxs already written
		std::vector<int> split_idxs_written;

		// Iterate over optimization steps
		for (int i_opt=0; i_opt<n_opt; i_opt++)
		{
			std::cout << "Opt step " << i_opt << " / " << n_opt-1 << std::endl;

			// Offset 
			i_opt_translated = i_opt + options.opt_idx_start_writing;

			// Exp decay
			if (options.exp_decay) {
				exp_decay_t0 = options.exp_decay_t0_values[i_opt];
				exp_decay_lambda = options.exp_decay_lambda_values[i_opt];
			};

			// Start and stop time
			if (options.time_cutoff) {
				t_start = options.time_cutoff_start_values[i_opt];
				t_end = options.time_cutoff_end_values[i_opt];
			} else {
				t_start = 0;
				t_end = _time.n();
			};

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
						itbf->nesterov_move_to_intermediate_pt(i_opt_translated);
					};
				};
			};

			// Write the basis funcs
			if (options.write && !options.write_bf_only_final) {
				write_bfs(options.dir_write+"F/",i_opt_translated);
			};

			/*****
			Step 2 - Iterate over batch
			*****/

			if (DIAG_SOLVE) { std::cout << "   Looping over batch" << std::endl; };

			// Reset writing
			split_idxs_written.clear();

			for (int i_batch=0; i_batch<fnames_to_use[i_opt].size(); i_batch++) 
			{
				if (options.verbose) {
					std::cout << "Doing sample: " << i_batch << " / " << fnames_to_use[i_opt].size() << " file: " << fnames_to_use[i_opt][i_batch].fname << std::endl;
				};

				/*****
				Step 2.1 - Read the IC
				*****/

				read_init_cond(fnames_to_use[i_opt][i_batch].fname_ic);

				/*****
				Step 2.2 - Solve the current trajectory
				*****/

				if (DIAG_SOLVE) { std::cout << "Solving ixn param" << std::endl; };

				solve_ixn_param_traj(t_end);

				// Write
				if (options.write && fnames_to_use[i_opt][i_batch].write && options.write_ixn_params) {
					auto f = std::find(split_idxs_written.begin(),split_idxs_written.end(),fnames_to_use[i_opt][i_batch].idx_split);
					if (f == split_idxs_written.end()) {
						write_ixn_params(options.dir_write+"ixn_params/",i_opt_translated,{fnames_to_use[i_opt][i_batch].idx_split});
						split_idxs_written.push_back(fnames_to_use[i_opt][i_batch].idx_split);
					};
				};

				if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

				/*****
				Step 2.3 - Solve the variational problem traj
				*****/

				if (DIAG_SOLVE) { std::cout << "Solving var term" << std::endl; };

				solve_var_traj();

				// Write
				if (options.write && options.write_var_terms && fnames_to_use[i_opt][i_batch].write) {
					write_var_terms(options.dir_write+"var_terms/",i_opt_translated);
				};
				
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
				for (_t_opt=t_start; _t_opt < _n_t_soln; _t_opt++)
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
						_latt.read_from_file(fnames_to_use[i_opt][i_batch].fname + pad_str(options.time_idx_start_reading+_t_opt,4) + ".txt");
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
								_latt.sample(_t_opt);
							} else {
								std::cerr << "Error! Probabilistic asleep visible units not supported yet!" << std::endl;
								exit(EXIT_FAILURE);
							};
						} else {
							if (options.asleep_final_visible_are_binary) {
								_latt.sample(_t_opt);
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

				if (options.write && fnames_to_use[i_opt][i_batch].write && options.write_moments) {
					write_moments(options.dir_write+"moments/",i_opt_translated,{fnames_to_use[i_opt][i_batch].idx_split,fnames_to_use[i_opt][i_batch].idx_sample});
				};

				/*****
				Step 2.7 - Gather the update (but dont commit)
				*****/

				for (auto itbf=_bfs.begin(); itbf!=_bfs.end(); itbf++) {
					itbf->update_gather(t_start, _n_t_soln, _time.delta(), dopt, options.exp_decay, exp_decay_t0, exp_decay_lambda, options.l2_reg_params_mode, l2_lambda_params, l2_reg_centers);
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
			write_bfs(options.dir_write+"F/",n_opt+options.opt_idx_start_writing);
		};
	};










	/****************************************
	Solve single IC, dividing up randomly the traj and restarting the diff eq integration for the var terms
	****************************************/

	void OptProblem::Impl::solve_rand_div(std::vector<std::string> fnames, int n_opt, int batch_size, int n_cd_steps, double dopt, int n_divisions, OptionsSolve options=OptionsSolve()) {

		// Pick filenames
		std::vector<std::vector<std::string>> fnames_to_use;
		std::vector<std::string> fnames_batch;
		std::vector<std::string> fnames_remaining;
		std::vector<std::string>::iterator itf;

		for (int i_opt=0; i_opt < n_opt; i_opt++) {

			// Clear batch
			fnames_batch.clear();

			// Reset remaining
			fnames_remaining=fnames;

			// Choose batch
			if (options.use_same_lattice_in_batch) {
				// Use a single lattice - pick a rand
				itf = fnames_remaining.begin();
				std::advance(itf,randI(0,fnames_remaining.size()-1));
				// Add
				fnames_batch.push_back(*itf);
				// Only this one
				for (int i_batch=1; i_batch<batch_size; i_batch++) {
					fnames_batch.push_back(fnames_batch.back());
				};
			} else {
				while (fnames_batch.size() < batch_size) {
					// Choose a random one
					itf = fnames_remaining.begin();
					std::advance(itf,randI(0,fnames_remaining.size()-1));
					// Add
					fnames_batch.push_back(*itf);
					// Don't choose again
					fnames_remaining.erase(itf);
				};
			};

			/*
			for (auto s: fnames_batch) {
				std::cout << s << " ";
			};
			std::cout << "" << std::endl;
			*/
			
			// Add the batch
			fnames_to_use.push_back(fnames_batch);
		};

		// Main solve function
		_solve_rand_div(fnames_to_use, n_opt, n_cd_steps, dopt, n_divisions, options);

	};
	void OptProblem::Impl::solve_rand_div(std::vector<std::vector<std::string>> fname_collection, int n_opt, int n_cd_steps, double dopt, int n_divisions, OptionsSolve options=OptionsSolve()) {

		// Pick filenames
		std::vector<std::vector<std::string>> fnames_to_use;
		std::vector<std::string> fnames_batch;
		std::vector<std::string>::iterator itf;

		for (int i_opt=0; i_opt < n_opt; i_opt++) {
			// Pick a random one
			fnames_batch = fname_collection[randI(0,fname_collection.size()-1)];

			// Add
			fnames_to_use.push_back(fnames_batch);
		};

		// Main solve function
		_solve_rand_div(fnames_to_use, n_opt, n_cd_steps, dopt, n_divisions, options);

	};
	void OptProblem::Impl::_solve_rand_div(std::vector<std::vector<std::string>> fnames_to_use, int n_opt, int n_cd_steps, double dopt, int n_divisions, OptionsSolve options=OptionsSolve()) {

		// Clear/make directories if needed
		if (options.clear_dir) {
			system(("rm -rf " + options.dir_write).c_str());
		};
		if (options.write) {
			system(("mkdir -p " + options.dir_write).c_str());
			system(("mkdir -p " + options.dir_write + "F").c_str());
			system(("mkdir -p " + options.dir_write + "ixn_params").c_str());
			system(("mkdir -p " + options.dir_write + "moments").c_str());
		};
		if (options.write_var_terms) {
			system(("mkdir -p " + options.dir_write + "var_terms").c_str());
		};

		// Write the grids
		if (options.write) {
			write_bf_grids(options.dir_write);
			write_t_grid(options.dir_write);
		};

		// Opt step with offset
		int i_opt;

		// Check exp decay / cutoff times
		if (options.exp_decay) {
			// Check
			if (options.exp_decay_t0_values.size() != n_opt || options.exp_decay_lambda_values.size() != n_opt) {
				std::cerr << "Error! In exp decay mode, must provide t0 and lambda values for all optimization timesteps." << std::endl;
				exit(EXIT_FAILURE);
			};
		};
		if (options.time_cutoff) {
			// Check
			if (options.time_cutoff_start_values.size() != n_opt || options.time_cutoff_end_values.size() != n_opt) {
				std::cerr << "Error! In time cutoff mode, must provide values for all optimization timesteps." << std::endl;
				exit(EXIT_FAILURE);
			};
		};

		// Convert l2 reg for params
		std::map<IxnParamTraj*,double> l2_lambda_params, l2_reg_centers;
		if (options.l2_reg_params_mode) {
			for (auto const &pr: options.l2_lambda_params) {
				l2_lambda_params[_not_nullptr(_find_ixn_param(pr.first))] = pr.second;
			};
			for (auto const &pr: options.l2_reg_centers) {
				l2_reg_centers[_not_nullptr(_find_ixn_param(pr.first))] = pr.second;
			};
		};

		// Values for exp decay
		double exp_decay_t0=0., exp_decay_lambda=0.;

		// Start/stop time for integration
		int t_start,t_end;

		// Check that no divisions specified is acceptable
		if (n_divisions > _time.n()) {
			std::cerr << "Error! specified: " << n_divisions << " divisions but there are only: " << _time.n() << " timesteps...." << std::endl;
			exit(EXIT_FAILURE);
		};

		// Struct to hold the chosen division timepoints
		std::vector<int> timepoints_all;
		for(auto i=0; i<_time.n(); i++) {
			timepoints_all.push_back(i);	
		};
		std::vector<int> division_times;

		// Iterate over optimization steps
		for (int i_opt_from_zero=0; i_opt_from_zero<n_opt; i_opt_from_zero++)
		{
			std::cout << "Opt step " << i_opt_from_zero << " / " << n_opt-1 << std::endl;

			// Offset
			i_opt = i_opt_from_zero + options.opt_idx_start_writing;

			// Exp decay
			if (options.exp_decay) {
				exp_decay_t0 = options.exp_decay_t0_values[i_opt_from_zero];
				exp_decay_lambda = options.exp_decay_lambda_values[i_opt_from_zero];
			};

			// Start and stop time
			if (options.time_cutoff) {
				t_start = options.time_cutoff_start_values[i_opt_from_zero];
				t_end = options.time_cutoff_end_values[i_opt_from_zero];
			} else {
				t_start = 0;
				t_end = _time.n();
			};

			/*****
			Step 0 - Check nesterov
			*****/

			if (options.nesterov) {
				// If first opt step, do nothing, but set the "prev" point to the current to initialize
				if (i_opt_from_zero == 0) {
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

			solve_ixn_param_traj(t_end);

			// Write
			if (options.write && options.write_ixn_params) {
				write_ixn_params(options.dir_write+"ixn_params/",i_opt);
			};

			/*****
			Step 2 - choose divisions
			*****/

			// Choose divisions
			std::random_shuffle ( timepoints_all.begin(), timepoints_all.end() );
			division_times = std::vector<int>(timepoints_all.begin(),timepoints_all.begin()+n_divisions-1);
			std::sort(division_times.begin(), division_times.end());
			// Adjust first point
			while (division_times.front() <= t_start) {
				division_times.erase(division_times.begin());
			};
			// Now the first point is greater than t_start
			division_times.insert(division_times.begin(),t_start);
			// Now the first point is t_start
			// Adjust last point
			while (division_times.back() >= _n_t_soln) {
				division_times.pop_back();
			};
			// Now the last point is less than _n_t_soln
			division_times.push_back(_n_t_soln);
			// Now the last point is _n_t_soln

			if (options.verbose) {
				std::cout << "Division times: " << std::flush;
				for (auto i=0; i<division_times.size(); i++) {
					std::cout << division_times[i] << " ";
				};
				std::cout << std::endl;
			};

			/*****
			Step 4 - reset the moments at all times
			*****/

			for (auto itp = _ixn_params.begin(); itp != _ixn_params.end(); itp++) {
				itp->moments_reset();
			};

			/*****
			Step 5 - Go through all times
			*****/

			// Current division
			int div_idx = -1;
			int t_div_start, t_div_end;
			// Start of the next div
			int t_div_start_next = t_start;

			// Use the class variable _t_opt to iterate
			for (_t_opt=t_start; _t_opt < _n_t_soln; _t_opt++)
			{
				if (options.verbose) {
					std::cout << "time: " << _t_opt << std::flush;
				};

				/*****
				Step 5.0 - division and var term
				*****/

				// Have we entered a new division?
				if (_t_opt == t_div_start_next) {

					// Yes, starting new div
					if (DIAG_SOLVE_RAND_DIV) {
						std::cout << "Starting new div..." << std::endl;
					};

					// First gather old update
					if (_t_opt != t_start) {
						if (DIAG_SOLVE_RAND_DIV) {
							std::cout << "Gathering update in div idx " << div_idx << " update from " << t_div_start << " to: " << t_div_end << std::endl;
						};
						for (auto itbf=_bfs.begin(); itbf!=_bfs.end(); itbf++) {
							itbf->update_gather(t_div_start, t_div_end, _time.delta(), dopt, options.exp_decay, exp_decay_t0, exp_decay_lambda, options.l2_reg_params_mode, l2_lambda_params, l2_reg_centers);
						};
					};

					// Advance the idxs
					div_idx++;
					t_div_start = t_div_start_next;
					t_div_end = division_times[div_idx+1];
					t_div_start_next = t_div_end;

					// Solve var trajs between these
					if (DIAG_SOLVE_RAND_DIV) { 
						std::cout << "Solving var term in div idx: " << div_idx << " times: " << t_div_start << " to: " << t_div_end << std::endl; 
					};

					solve_var_traj_from_zero(t_div_start,t_div_end);

					// Write
					/*
					if (options.write && options.write_var_terms) {
						write_var_terms(options.dir_write+"var_terms/",i_opt);
					};
					*/
				};

				if (options.verbose) {
					std::cout << " div: " << div_idx << std::flush;
				};

				/*****
				Step 5.1 - loop over all samples in the batch
				*****/

				for (int i_batch=0; i_batch<fnames_to_use[i_opt_from_zero].size(); i_batch++) 
				{
					if (options.verbose) {
						std::cout << "." << std::flush;
					};

					/*****
					Step 5.1.1 - Read in sample at this timestep
					*****/

					if (options.awake_visible_are_binary) {
						// Binary
						_latt.read_from_file(fnames_to_use[i_opt_from_zero][i_batch] + pad_str(options.time_idx_start_reading+_t_opt,4) + ".txt");
					} else {
						// Probabilistic
						std::cerr << "Error! Probabilistic awake visible units not supported yet!" << std::endl;
						exit(EXIT_FAILURE);
					};
 
					/*****
					Step 5.1.2 - Activate the hidden units if needed
					*****/

					if (_hidden_layer_exists) {
						for (auto ithu = _hidden_units.begin(); ithu != _hidden_units.end(); ithu++) {
							ithu->activate(options.awake_hidden_are_binary);
						};
					};

					/*****
					Step 5.1.3 - Record the awake moments
					*****/

					for (auto itp = _ixn_params.begin(); itp != _ixn_params.end(); itp++) {
						itp->moments_retrieve_at_time(IxnParamTraj::AWAKE,_t_opt,fnames_to_use[i_opt_from_zero].size());
					};

					/*****
					Step 5.1.4 - CD steps
					*****/

					for (int cd_step=0; cd_step<n_cd_steps; cd_step++)
					{
						// Sample

						if (cd_step != n_cd_steps-1) {
							if (options.asleep_visible_are_binary) {
								_latt.sample(_t_opt);
							} else {
								std::cerr << "Error! Probabilistic asleep visible units not supported yet!" << std::endl;
								exit(EXIT_FAILURE);
							};
						} else {
							if (options.asleep_final_visible_are_binary) {
								_latt.sample(_t_opt);
							} else {
								std::cerr << "Error! Probabilistic asleep visible units not supported yet!" << std::endl;
								exit(EXIT_FAILURE);
							};
						};

						// Activate the hidden units if needed

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

					for (auto itp = _ixn_params.begin(); itp != _ixn_params.end(); itp++) {
						itp->moments_retrieve_at_time(IxnParamTraj::ASLEEP,_t_opt,fnames_to_use[i_opt_from_zero].size());
					};
				};

				// If the last timestep, also gather the update
				if (_t_opt == _n_t_soln-1) {
					if (DIAG_SOLVE_RAND_DIV) {
						std::cout << "Gathering update in div idx " << div_idx << " update from " << t_div_start << " to: " << t_div_end << std::endl;
					};
					for (auto itbf=_bfs.begin(); itbf!=_bfs.end(); itbf++) {
						itbf->update_gather(t_div_start, t_div_end, _time.delta(), dopt, options.exp_decay, exp_decay_t0, exp_decay_lambda, options.l2_reg_params_mode, l2_lambda_params, l2_reg_centers);
					};
				};

				if (options.verbose) {
					std::cout << std::endl;
				};
			};

			/*****
			Step 6 - Commit gathered update
			*****/

			for (auto itbf=_bfs.begin(); itbf!=_bfs.end(); itbf++) {
				itbf->update_committ_gathered();
			};

			/*****
			Step 7 - Write if needed
			*****/

			// Write the moments if needed
			if (options.write && options.write_moments) {
				write_moments(options.dir_write+"moments/",i_opt);
			};

			// Write var terms if needed
			if (options.write && options.write_var_terms) {
				write_var_terms(options.dir_write+"var_terms/",i_opt);
			};
		};

		// Write the basis funcs one last time
		if (options.write) {
			write_bfs(options.dir_write+"F/",n_opt+options.opt_idx_start_writing);
		};
		
	};


















	/********************
	Read some initial conditions
	********************/

	void OptProblem::Impl::read_init_cond(std::string fname_ic) {

		std::ifstream f;
		f.open(fname_ic, std::ios::in);
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

	void OptProblem::Impl::write_ixn_params(std::string dir, int idx_opt_step) const {
		for (auto it = _ixn_params.begin(); it!=_ixn_params.end(); it++) {
			it->write_vals(dir,idx_opt_step,{},_n_t_soln);
		};
	};
	void OptProblem::Impl::write_ixn_params(std::string dir, int idx_opt_step, std::vector<int> idxs) const {
		for (auto it = _ixn_params.begin(); it!=_ixn_params.end(); it++) {
			it->write_vals(dir,idx_opt_step,idxs,_n_t_soln);
		};
	};
	void OptProblem::Impl::write_bfs(std::string dir, int idx_opt_step) const {
		for (auto it=_bfs.begin(); it!=_bfs.end(); it++) {
			it->write_vals(dir, idx_opt_step);
		};
	};
	void OptProblem::Impl::write_var_terms(std::string dir, int idx_opt_step) const {
		for (auto it=_var_terms.begin(); it!=_var_terms.end(); it++) {
			it->write_vals(dir, idx_opt_step);
		};
	};
	void OptProblem::Impl::write_moments(std::string dir, int idx_opt_step) const {
		for (auto it = _ixn_params.begin(); it!=_ixn_params.end(); it++) {
			it->write_moments(dir,idx_opt_step,{},_n_t_soln);
		};
	};
	void OptProblem::Impl::write_moments(std::string dir, int idx_opt_step, std::vector<int> idxs) const {
		for (auto it = _ixn_params.begin(); it!=_ixn_params.end(); it++) {
			it->write_moments(dir,idx_opt_step,idxs,_n_t_soln);
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

	void OptProblem::Impl::read_init_cond_for_ixn_param(std::string ixn_func_name, std::string fname, int line_idx) {
		// Find the ixn func
		IxnParamTraj* ixn_func = _find_ixn_param(ixn_func_name);
		if (!ixn_func) {
			std::cerr << "ERROR: Could not find ixn func to read init cond into." << std::endl;
			exit(EXIT_FAILURE);
		};
		// Read
		std::ifstream f;
		f.open(fname);
		std::string line;
		std::string s="";
		std::istringstream iss;
		double init_cond;
		int i_line = 0;
		if (f.is_open()) { // make sure we found it
			while (getline(f,line)) {
				if (line == "") { continue; };
				iss = std::istringstream(line);
				iss >> s;
				init_cond = atof(s.c_str());
				s=""; // reset
				// Check
				if (i_line == line_idx) {
					// Set latest
					ixn_func->set_init_cond(init_cond);
					// Fin
					return;
				};
				// Next line
				i_line++;
			};
		};
		// Set latest
		ixn_func->set_init_cond(init_cond);
	};

	/****************************************
	OptProblem - IMPLEMENTATION - PRIVATE
	****************************************/

	/********************
	Search functions
	********************/

	Species* OptProblem::Impl::_not_nullptr(Species* ptr) {
		if (!ptr) {
			std::cerr << "ERROR: could not find species!" << std::endl;
			exit(EXIT_FAILURE);
		};
		return ptr;
	};
	HiddenSpecies* OptProblem::Impl::_not_nullptr(HiddenSpecies* ptr) {
		if (!ptr) {
			std::cerr << "ERROR: could not find hidden species!" << std::endl;
			exit(EXIT_FAILURE);
		};
		return ptr;
	};
	IxnParamTraj* OptProblem::Impl::_not_nullptr(IxnParamTraj* ptr) {
		if (!ptr) {
			std::cerr << "ERROR: could not find IxnParamTraj!" << std::endl;
			exit(EXIT_FAILURE);
		};
		return ptr;
	};
	BasisFunc* OptProblem::Impl::_not_nullptr(BasisFunc* ptr) {
		if (!ptr) {
			std::cerr << "ERROR: could not find BasisFunc!" << std::endl;
			exit(EXIT_FAILURE);
		};
		return ptr;
	};
	VarTermTraj* OptProblem::Impl::_not_nullptr(VarTermTraj* ptr) {
		if (!ptr) {
			std::cerr << "ERROR: could not find VarTermTraj!" << std::endl;
			exit(EXIT_FAILURE);
		};
		return ptr;
	};
	HiddenSpecies* OptProblem::Impl::_find_hidden_species(std::string name) {
		for (auto it=_hidden_species.begin(); it!=_hidden_species.end(); it++) {
			if (it->name() == name) {
				return &*it;
			};
		};
		return nullptr;
	};
	Species* OptProblem::Impl::_find_species(std::string name) {
		for (auto it=_species.begin(); it!=_species.end(); it++) {
			if (it->name() == name) {
				return &*it;
			};
		};
		return nullptr;
	};
	IxnParamTraj* OptProblem::Impl::_find_ixn_param(std::string name) {
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			if (it->name() == name) {
				return &*it;
			};
		};
		return nullptr;
	};
	BasisFunc* OptProblem::Impl::_find_basis_func(std::string name) {
		for (auto it=_bfs.begin(); it!=_bfs.end(); it++) {
			if (it->name() == name) {
				return &*it;
			};
		};
		return nullptr;
	};
	VarTermTraj* OptProblem::Impl::_find_var_term(std::string name) {
		for (auto it=_var_terms.begin(); it!=_var_terms.end(); it++) {
			if (it->name() == name) {
				return &*it;
			};
		};
		return nullptr;
	};
	Counter* OptProblem::Impl::_find_ctr_by_species(std::string s) {
		for (auto it=_counters.begin(); it!=_counters.end(); it++) {
			if (it->is_counting_species(s)) {
				return &*it;
			};
		};
		return nullptr;
	};
	Counter* OptProblem::Impl::_find_ctr_by_species(std::string s1, std::string s2) {
		for (auto it=_counters.begin(); it!=_counters.end(); it++) {
			if (it->is_counting_species(s1,s2)) {
				return &*it;
			};
		};
		return nullptr;
	};
	Counter* OptProblem::Impl::_find_ctr_by_species(std::string s1, std::string s2, std::string s3) {
		for (auto it=_counters.begin(); it!=_counters.end(); it++) {
			if (it->is_counting_species(s1,s2,s3)) {
				return &*it;
			};
		};
		return nullptr;
	};
	Counter* OptProblem::Impl::_find_ctr_by_species(std::string s1, std::string s2, std::string s3, std::string s4) {
		for (auto it=_counters.begin(); it!=_counters.end(); it++) {
			if (it->is_counting_species(s1,s2,s3,s4)) {
				return &*it;
			};
		};
		return nullptr;
	};

	/********************
	Add counter
	********************/

	void OptProblem::Impl::_add_counter(std::string s) {
		// Check it does not exist
		Counter *ctr = _find_ctr_by_species(s);
		if (ctr) {
			return;
		};

		// Find species
		Species* sp = _not_nullptr(_find_species(s));

		// Make counter
		_counters.push_back(Counter(sp));
		
		// Add this counter to the species
		sp->set_counter(&_counters.back());
	};
	void OptProblem::Impl::_add_counter(std::string s1, std::string s2) {
		// Check it does not exist
		Counter *ctr = _find_ctr_by_species(s1,s2);
		if (ctr) {
			return;
		};

		// Find species
		Species* sp1 = _not_nullptr(_find_species(s1));
		Species* sp2 = _not_nullptr(_find_species(s2));

		// Make counter
		_counters.push_back(Counter(sp1,sp2));
		
		// Add this counter to the species
		sp1->add_nn_counter(sp2,&_counters.back());
		sp2->add_nn_counter(sp1,&_counters.back());
	};
	void OptProblem::Impl::_add_counter(std::string s1, std::string s2, std::string s3) {
		// Check it does not exist
		Counter *ctr = _find_ctr_by_species(s1,s2,s3);
		if (ctr) {
			return;
		};

		// Find species
		Species* sp1 = _not_nullptr(_find_species(s1));
		Species* sp2 = _not_nullptr(_find_species(s2));
		Species* sp3 = _not_nullptr(_find_species(s3));

		// Make counter
		_counters.push_back(Counter(sp1,sp2,sp3));
		
		// Add this counter to the species
		sp1->add_triplet_counter(sp2,sp3,&_counters.back());
		sp2->add_triplet_counter(sp1,sp3,&_counters.back());
		sp3->add_triplet_counter(sp1,sp2,&_counters.back());
	};
	void OptProblem::Impl::_add_counter(std::string s1, std::string s2, std::string s3, std::string s4) {
		// Check it does not exist
		Counter *ctr = _find_ctr_by_species(s1,s2,s3,s4);
		if (ctr) {
			return;
		};

		// Find species
		Species* sp1 = _not_nullptr(_find_species(s1));
		Species* sp2 = _not_nullptr(_find_species(s2));
		Species* sp3 = _not_nullptr(_find_species(s3));
		Species* sp4 = _not_nullptr(_find_species(s4));

		// Make counter
		_counters.push_back(Counter(sp1,sp2,sp3,sp4));
		
		// Add this counter to the species
		sp1->add_quartic_counter(sp2,sp3,sp4,&_counters.back());
		sp2->add_quartic_counter(sp1,sp3,sp4,&_counters.back());
		sp3->add_quartic_counter(sp1,sp2,sp4,&_counters.back());
		sp4->add_quartic_counter(sp1,sp2,sp3,&_counters.back());
	};










































	/****************************************
	OptProblem IMPL forwards
	****************************************/

	// Constructor
	OptProblem::OptProblem(std::vector<Dim> dims, std::vector<std::string> species_visible, std::vector<std::string> species_hidden, double t_max, int n_t, int box_length, int lattice_dim) : _impl(new Impl(dims,species_visible,species_hidden,t_max,n_t,box_length,lattice_dim)) {};
	OptProblem::OptProblem(const OptProblem& other) : _impl(new Impl(*other._impl)) {};
	OptProblem::OptProblem(OptProblem&& other) = default;
	OptProblem& OptProblem::operator=(OptProblem other) {
        _impl = std::move(other._impl);
        return *this; 
	};
	OptProblem::~OptProblem() = default;

	/********************
	Change the time limit
	********************/

	void OptProblem::set_n_t(int n_t) {
		_impl->set_n_t(n_t);
	};

	/********************
	Set IC for ixn param
	********************/

	void OptProblem::set_ic_for_ixn_param(std::string param_name, double val) {
		_impl->set_ic_for_ixn_param(param_name,val);
	};

	void OptProblem::set_fixed_awake_moment_for_dim(std::string param_name, std::vector<double> vals) {
		_impl->set_fixed_awake_moment_for_dim(param_name,vals);
	};

	void OptProblem::add_hidden_unit(std::vector<std::string> species_possible, std::vector<std::vector<int>> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params) {
		_impl->add_hidden_unit(species_possible,lattice_idxs,w_params,b_params);
	};
	void OptProblem::add_hidden_unit(std::vector<std::string> species_possible, std::vector<int> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params) {
		_impl->add_hidden_unit(species_possible,lattice_idxs,w_params,b_params);
	};

	void OptProblem::validate_setup() const {
		_impl->validate_setup();
	};
	void OptProblem::validate_graph() const {
		_impl->validate_graph();
	};

	void OptProblem::solve_ixn_param_traj() {
		_impl->solve_ixn_param_traj();
	};
	void OptProblem::solve_var_traj() {
		_impl->solve_var_traj();
	};
	void OptProblem::solve_var_traj_from_zero(int it_start, int it_end) {
		_impl->solve_var_traj_from_zero(it_start,it_end);
	};

	void OptProblem::solve(std::vector<std::string> fnames, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolve options) {
		_impl->solve(fnames,n_opt,batch_size,n_cd_steps,dopt,options);
	};
	void OptProblem::solve(std::vector<std::vector<std::string>> fname_collection, int n_opt, int n_cd_steps, double dopt, OptionsSolve options) {
		_impl->solve(fname_collection,n_opt,n_cd_steps,dopt,options);
	};

	void OptProblem::solve_varying_ic(std::vector<FName> fnames, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolve options) {
		_impl->solve_varying_ic(fnames,{},n_opt,batch_size,n_cd_steps,dopt,options);
	};
	void OptProblem::solve_varying_ic(std::vector<FName> fnames, std::vector<FName> fnames_used_in_every_batch, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolve options) {
		_impl->solve_varying_ic(fnames,fnames_used_in_every_batch,n_opt,batch_size,n_cd_steps,dopt,options);
	};
	void OptProblem::solve_varying_ic(std::vector<std::vector<FName>> fname_collection, int n_opt, int n_cd_steps, double dopt, OptionsSolve options) {
		_impl->solve_varying_ic(fname_collection,n_opt,n_cd_steps,dopt,options);
	};

	void OptProblem::solve_rand_div(std::vector<std::string> fnames, int n_opt, int batch_size, int n_cd_steps, double dopt, int n_divisions, OptionsSolve options) {
		_impl->solve_rand_div(fnames,n_opt,batch_size,n_cd_steps,dopt,n_divisions,options);
	};
	void OptProblem::solve_rand_div(std::vector<std::vector<std::string>> fnames, int n_opt, int n_cd_steps, double dopt, int n_divisions, OptionsSolve options) {
		_impl->solve_rand_div(fnames,n_opt,n_cd_steps,dopt,n_divisions,options);
	};

	void OptProblem::read_basis_func(std::string bf_name, std::string fname) {
		_impl->read_basis_func(bf_name,fname);
	};

	void OptProblem::read_init_cond_for_ixn_param(std::string ixn_func_name, std::string fname, int line_idx) {
		_impl->read_init_cond_for_ixn_param(ixn_func_name, fname, line_idx);
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
		_impl->write_ixn_params(dir,idx1,std::vector<int>({idx2}));
	};
	void OptProblem::write_ixn_params(std::string dir, int idx, std::vector<int> idxs) const {
		_impl->write_ixn_params(dir,idx,idxs);
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
		_impl->write_moments(dir,idx1,std::vector<int>({idx2}));
	};
	void OptProblem::write_moments(std::string dir, int idx, std::vector<int> idxs) const {
		_impl->write_moments(dir,idx,idxs);
	};

};





