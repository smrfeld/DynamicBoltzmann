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
	Dim - IMPLEMENTATION
	****************************************/

	class Dim::Impl {
	private:

		// Name
		std::string _name;

		// Type
		DimType _type;

		// Name of associated species
		bool _any_species;
		std::vector<std::string> _species;
		std::vector<std::vector<std::string> > _species_multiple;

		// Min/max/npts/initial value
		double _min,_max;
		int _n;
		double _init;

		// Basis function dimension names
		std::vector<std::string> _basis_func_dims;

		// Constructor helpers
		void _clean_up();
		void _copy(const Impl& other);
		void _reset();
		void _shared_constructor(std::string name, DimType type, std::vector<std::string> basis_func_dims, double min, double max, int n, double init);

	public:

		// Constructor
		Impl(std::string name, DimType type, std::vector<std::string> basis_func_dims, double min, double max, int n, double init);
		Impl(std::string name, DimType type, std::string species, std::vector<std::string> basis_func_dims, double min, double max, int n, double init);
		Impl(std::string name, DimType type, std::vector<std::string> species, std::vector<std::string> basis_func_dims, double min, double max, int n, double init);
		Impl(std::string name, DimType type, std::vector<std::vector<std::string>> species, std::vector<std::string> basis_func_dims, double min, double max, int n, double init);
		Impl(const Impl& other);
		Impl(Impl&& other);
	    Impl& operator=(const Impl& other);
	    Impl& operator=(Impl&& other);
		~Impl();

		/********************
		Getters
		********************/

		// Name
		std::string name() const;

		// Type
		DimType type() const;

		// Does it apply to all species?
		bool any_species() const;

		// Basis func dims
		std::vector<std::string> basis_func_dims() const;

		// Min/max/n/init
		double min() const;
		double max() const;
		double n() const;
		double init() const;

		// Get species
		std::vector<std::string> get_species_h() const;
		std::vector<std::string> get_species_b() const;
		std::vector<std::vector<std::string>> get_species_J() const;
		std::vector<std::vector<std::string>> get_species_K() const;
		std::vector<std::vector<std::string>> get_species_W() const;

		/********************
		Setters
		********************/

		// Add basis func dimension
		void add_basis_func_dim(std::string basis_func);

		// Add species
		void add_species_h(std::string species);
		void add_species_b(std::string species);
		void add_species_J(std::string species1, std::string species2);
		void add_species_K(std::string species1, std::string species2, std::string species3);
		void add_species_W(std::string species_visible, std::string species_hidden);
	};





































	/****************************************
	Dim - IMPLEMENTATION DEFINITIONS
	****************************************/

	/********************
	Constructor
	********************/

	Dim::Impl::Impl(std::string name, DimType type, std::vector<std::string> basis_func_dims, double min, double max, int n, double init) {
		_shared_constructor(name, type, basis_func_dims, min, max, n, init);
		// Yes any species
		_any_species = true;
	};
	Dim::Impl::Impl(std::string name, DimType type, std::string species, std::vector<std::string> basis_func_dims, double min, double max, int n, double init) {
		_shared_constructor(name, type, basis_func_dims, min, max, n, init);
		// Not any species
		_any_species = false;
		// Add
		if (_type == H) {
			add_species_h(species);
		} else if (_type == B) {
			add_species_b(species);
		};
	};
	Dim::Impl::Impl(std::string name, DimType type, std::vector<std::string> species, std::vector<std::string> basis_func_dims, double min, double max, int n, double init) {
		_shared_constructor(name, type, basis_func_dims, min, max, n, init);
		// Not any species
		_any_species = false;
		// Add
		if (_type == H) {
			for (auto s: species) {
				add_species_h(s);
			};
		} else if (_type == B) {
			for (auto s: species) {
				add_species_b(s);
			};
		} else if (type == J) {
			if (species.size() != 2) {
				std::cerr << "Error! must be 2 species for J" << std::endl;
				exit(EXIT_FAILURE);
			};
			add_species_J(species[0],species[1]);
		} else if (type == W) {
			if (species.size() != 2) {
				std::cerr << "Error! must be 2 species for W" << std::endl;
				exit(EXIT_FAILURE);
			};
			add_species_W(species[0],species[1]);
		};
	};
	Dim::Impl::Impl(std::string name, DimType type, std::vector<std::vector<std::string>> species, std::vector<std::string> basis_func_dims, double min, double max, int n, double init) {
		// Check type
		if (type == B || type != H) {
			std::cerr << "Error! Only for J or W." << std::endl;
			exit(EXIT_FAILURE);
		};
		_shared_constructor(name, type, basis_func_dims, min, max, n, init);
		// Not any species
		_any_species = false;
		// Add
		if (type == J) {
			for (auto s_pair: species) {
				if (s_pair.size() != 2) {
					std::cerr << "Error! must be 2 species for J" << std::endl;
					exit(EXIT_FAILURE);
				};
				add_species_J(s_pair[0],s_pair[1]);
			};
		} else if (type == W) {
			for (auto s_pair: species) {
				if (s_pair.size() != 2) {
					std::cerr << "Error! must be 2 species for W" << std::endl;
					exit(EXIT_FAILURE);
				};
				add_species_W(s_pair[0],s_pair[1]);
			};
		};
	};
	Dim::Impl::Impl(const Impl& other) {
		_copy(other);
	};
	Dim::Impl::Impl(Impl&& other) {
		_copy(other);
		other._reset();
	};
    Dim::Impl& Dim::Impl::operator=(const Impl& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    Dim::Impl& Dim::Impl::operator=(Impl&& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
    };
	Dim::Impl::~Impl()
	{
		_clean_up();
	};
	void Dim::Impl::_clean_up() {
		// Nothing...
	};
	void Dim::Impl::_copy(const Impl& other) {
		_shared_constructor(other._name, other._type,other._basis_func_dims, other._min, other._max, other._n, other._init);
		_any_species = other._any_species;
		_species = other._species;
		_species_multiple = other._species_multiple;
	};
	void Dim::Impl::_reset() {
		_name = "";
		_basis_func_dims.clear();
		_any_species = false;
		_species.clear();
		_species_multiple.clear();
		_min = 0.;
		_max = 0.;
		_n = 0;
		_init = 0.;
	};
	void Dim::Impl::_shared_constructor(std::string name, DimType type, std::vector<std::string> basis_func_dims, double min, double max, int n, double init) {
		_name = name;
		_type = type;
		_basis_func_dims = basis_func_dims;
		_min = min;
		_max = max;
		_n = n;
		_init = init;
	};

	/********************
	Getters
	********************/

	// Name
	std::string Dim::Impl::name() const {
		return _name;
	};

	// Type
	DimType Dim::Impl::type() const {
		return _type;
	};

	// Does it apply to all species?
	bool Dim::Impl::any_species() const {
		return _any_species;
	};

	// Basis func dims
	std::vector<std::string> Dim::Impl::basis_func_dims() const {
		return _basis_func_dims;
	};

	// Min/max/n/init
	double Dim::Impl::min() const {
		return _min;
	};
	double Dim::Impl::max() const {
		return _max;
	};
	double Dim::Impl::n() const {
		return _n;
	};
	double Dim::Impl::init() const {
		return _init;
	};

	// Get species
	std::vector<std::string> Dim::Impl::get_species_h() const {
		if (_type != H) {
			std::cerr << "Error! Requested species but not of type H." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _species;
	};
	std::vector<std::string> Dim::Impl::get_species_b() const {
		if (_type != B) {
			std::cerr << "Error! Requested species but not of type B." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _species;
	};
	std::vector<std::vector<std::string>> Dim::Impl::get_species_J() const {
		if (_type != J) {
			std::cerr << "Error! Requested species but not of type J." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _species_multiple;
	};
	std::vector<std::vector<std::string>> Dim::Impl::get_species_W() const {
		if (_type != W) {
			std::cerr << "Error! Requested species but not of type W." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _species_multiple;
	};

	/********************
	Setters
	********************/

	// Add basis func dimension
	void Dim::Impl::add_basis_func_dim(std::string dim) {
		_basis_func_dims.push_back(dim);
	};

	// Add species
	void Dim::Impl::add_species_h(std::string species) {
		// Check type
		if (_type != H) {
			std::cerr << "Error! Not H." << std::endl;
			exit(EXIT_FAILURE);
		};
		_any_species = false;
		_species.push_back(species);
	};
	void Dim::Impl::add_species_b(std::string species) {
		// Check type
		if (_type != B) {
			std::cerr << "Error! Not B." << std::endl;
			exit(EXIT_FAILURE);
		};
		_any_species = false;
		_species.push_back(species);
	};
	void Dim::Impl::add_species_J(std::string species1, std::string species2) {
		// Check type
		if (_type != J) {
			std::cerr << "Error! Not J." << std::endl;
			exit(EXIT_FAILURE);
		};
		_any_species = false;
		_species_multiple.push_back(std::vector<std::string>({species1,species2}));
	};
	void Dim::Impl::add_species_W(std::string species_visible, std::string species_hidden) {
		// Check type
		if (_type != W) {
			std::cerr << "Error! Not W." << std::endl;
			exit(EXIT_FAILURE);
		};
		_any_species = false;
		_species_multiple.push_back(std::vector<std::string>({species_visible,species_hidden}));
	};
































	/****************************************
	Dim - Impl forwards
	****************************************/

	Dim::Dim(std::string name, DimType type, std::vector<std::string> basis_func_dims, double min, double max, int n, double init) : _impl(new Impl(name,type,basis_func_dims,min,max,n,init)) {};
	Dim::Dim(std::string name, DimType type, std::string species, std::vector<std::string> basis_func_dims, double min, double max, int n, double init) : _impl(new Impl(name,type,species,basis_func_dims,min,max,n,init)) {};
	Dim::Dim(std::string name, DimType type, std::vector<std::string> species, std::vector<std::string> basis_func_dims, double min, double max, int n, double init) : _impl(new Impl(name,type,species,basis_func_dims,min,max,n,init)) {};
	Dim::Dim(std::string name, DimType type, std::vector<std::vector<std::string>> species, std::vector<std::string> basis_func_dims, double min, double max, int n, double init) : _impl(new Impl(name,type,species,basis_func_dims,min,max,n,init)) {};
	Dim::Dim(const Dim& other) : _impl(new Impl(*other._impl)) {};
	Dim::Dim(Dim&& other) = default;
	Dim& Dim::operator=(Dim other) {
        _impl = std::move(other._impl);
        return *this; 
	};
	Dim::~Dim() = default;

	// Name
	std::string Dim::name() const {
		return _impl->name();
	};

	// Type
	DimType Dim::type() const {
		return _impl->type();
	};

	// Does it apply to all species?
	bool Dim::any_species() const {
		return _impl->any_species();
	};

	// Basis func dims
	std::vector<std::string> Dim::basis_func_dims() const {
		return _impl->basis_func_dims();
	};

	// Min/max/n/init
	double Dim::min() const {
		return _impl->min();
	};
	double Dim::max() const {
		return _impl->max();
	};
	double Dim::n() const {
		return _impl->n();
	};
	double Dim::init() const {
		return _impl->init();
	};

	// Get species
	std::vector<std::string> Dim::get_species_h() const {
		return _impl->get_species_h();
	};
	std::vector<std::string> Dim::get_species_b() const {
		return _impl->get_species_b();
	};
	std::vector<std::vector<std::string>> Dim::get_species_J() const {
		return _impl->get_species_J();
	};
	std::vector<std::vector<std::string>> Dim::get_species_W() const {
		return _impl->get_species_W();
	};

	/********************
	Setters
	********************/

	// Add basis func dimension
	void Dim::add_basis_func_dim(std::string dim) {
		_impl->add_basis_func_dim(dim);
	};

	// Add species
	void Dim::add_species_h(std::string species) {
		_impl->add_species_h(species);
	};
	void Dim::add_species_b(std::string species) {
		_impl->add_species_b(species);
	};
	void Dim::add_species_J(std::string species1, std::string species2) {
		_impl->add_species_J(species1,species2);
	};
	void Dim::add_species_W(std::string species_visible, std::string species_hidden) {
		_impl->add_species_W(species_visible,species_hidden);
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

		Impl(std::vector<Dim> dims, std::vector<std::string> species_visible, std::vector<std::string> species_hidden, double t_max, int n_t, int box_length, int lattice_dim);
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
		void solve(std::vector<std::vector<std::string>> fnames, int n_opt, int n_cd_steps, double dopt, OptionsSolve options);
		void _solve(std::vector<std::vector<std::string>> fnames_to_use, int n_opt, int n_cd_steps, double dopt, OptionsSolve options);

		void solve_varying_ic(std::vector<FName> fnames, std::vector<FName> fnames_used_in_every_batch, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolve options);
		void solve_varying_ic(std::vector<std::vector<FName>> fnames, int n_opt, int n_cd_steps, double dopt, OptionsSolve options);
		void _solve_varying_ic(std::vector<std::vector<FName>> fnames_to_use, int n_opt, int n_cd_steps, double dopt, OptionsSolve options);

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

		// Tell the lattice about what dims exist
		if (DIAG_SETUP) { std::cout << "Telling lattice what dims exist..." << std::flush; };
		for (auto const &d: dims) {
			if (d.type()==H) {
				_latt.set_exists_h(true);
			} else if (d.type()==J) {
				_latt.set_exists_j(true);
			} else if (d.type()==W) {
				_hidden_layer_exists = true;
				_latt.set_exists_w(true);
			} else if (d.type()==B) {
				_hidden_layer_exists = true;
				// No need to tell lattice, since it only affects hidden unit activation and not sampling
			};
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Create the visible species
		if (DIAG_SETUP) { std::cout << "Create species..." << std::flush; };
		for (auto const &s: species_visible) {
			_species.push_back(Species(s));

			// Add to the lattice
			_latt.add_species(&(_species.back()));

			// Add the optimization time
			_species.back().set_opt_time_ptr(&_t_opt);
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Create the interaction params
		if (DIAG_SETUP) { std::cout << "Create ixn params..." << std::flush; };
		std::vector<std::string> s_names;
		std::vector<std::vector<std::string>> ss_names;
		for (auto const &d: dims) 
		{
			if (d.type()==H) {
				s_names = d.get_species_h();
				_ixn_params.push_back(IxnParamTraj(d.name(),Hp,_find_species(s_names[0]),d.min(),d.max(),d.n(),d.init(),n_t));
			} else if (d.type()==J) { 
				ss_names = d.get_species_J();
				_ixn_params.push_back(IxnParamTraj(d.name(),Jp,_find_species(ss_names[0][0]),_find_species(ss_names[0][1]),d.min(),d.max(),d.n(),d.init(),n_t));
			} else if (d.type()==W) { 
				ss_names = d.get_species_W();
				_ixn_params.push_back(IxnParamTraj(d.name(),Wp,_find_species(ss_names[0][0]),d.min(),d.max(),d.n(),d.init(),n_t));
			} else if (d.type()==B) {
				s_names = d.get_species_b();
				_ixn_params.push_back(IxnParamTraj(d.name(),Bp,_find_species(s_names[0]),d.min(),d.max(),d.n(),d.init(),n_t));
			};
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };

		// Add the interaction params to the species
		if (DIAG_SETUP) { std::cout << "Add ixn params to species..." << std::flush; };
		Species *sp=nullptr, *sp1=nullptr, *sp2=nullptr;
		IxnParamTraj *ip_ptr=nullptr;
		for (auto const &d: dims) {
			ip_ptr = _find_ixn_param_by_name(d.name());
			if (d.type()==H) {
				s_names = d.get_species_h();
				sp = _find_species(s_names[0]);
				sp->set_h_ptr(ip_ptr);
			} else if (d.type()==J) {
				ss_names = d.get_species_J();
				sp1 = _find_species(ss_names[0][0]);
				sp2 = _find_species(ss_names[0][1]);
				sp1->add_j_ptr(sp2,ip_ptr);
				sp2->add_j_ptr(sp1,ip_ptr);
			} else if (d.type()==W) {
				ss_names = d.get_species_W();
				sp = _find_species(ss_names[0][0]);
				sp->set_w_ptr(ip_ptr);		
			};
			// No need to tell species about biases
		};
		if (DIAG_SETUP) { std::cout << "ok." << std::endl; };


		// Ensure the J of the species are complete - if there is no ixn param that descripes the coupling, add a nullptr entry in the dictionary - later check if nullptr, then return 0
		// I think this is faster - otherwise there would be no reason to do it
		if (DIAG_SETUP) { std::cout << "Ensuring J are complete..." << std::flush; };
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
		for (auto const &d: dims) {
			// Find the basis func dimensions
			bf_ips.clear();
			for (auto bfd: d.basis_func_dims()) {
				bf_ips.push_back(_find_ixn_param_by_name(bfd));
			};
			// Make the basis function
			_bfs.push_back(BasisFunc("F_"+d.name(),bf_ips));
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
		for (auto const &num: dims) {
			// Go through denominators
			for (auto const &denom: dims) {
				// Find the basis func dimensions
				bf_ips.clear();
				for (auto bfd: denom.basis_func_dims()) {
					bf_ips.push_back(_find_ixn_param_by_name(bfd));
				};
				// Find the basis func
				bf_ptr = _find_basis_func("F_"+denom.name());
				// Find the interaction param
				ixn_param_ptr = _find_ixn_param_by_name(num.name());
				// Find the basis func corresponding to the numerator
				num_bf_ptr = _find_basis_func("F_"+num.name());
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
					bf_ips.push_back(_find_ixn_param_by_name(bfd));
				};
				// Find the variational term
				vt_ptr = _find_var_term("var_"+num.name()+"_wrt_"+denom->name());
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

		// Iterate over optimization steps
		for (int i_opt_from_zero=0; i_opt_from_zero<n_opt; i_opt_from_zero++)
		{
			std::cout << "Opt step " << i_opt_from_zero << " / " << n_opt-1 << std::endl;

			// Offset
			i_opt = i_opt_from_zero + options.opt_idx_start_writing;

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
			if (options.write_var_terms) {
				write_var_terms(options.dir_write+"var_terms/",i_opt);
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

				for (int i_batch=0; i_batch<fnames_to_use[_t_opt].size(); i_batch++) 
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
						_latt.read_from_file(fnames_to_use[_t_opt][i_batch] + pad_str(options.time_idx_start_reading+_t_opt,4) + ".txt");
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
						itp->moments_retrieve_at_time(IxnParamTraj::AWAKE,_t_opt,fnames_to_use[_t_opt].size());
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
						itp->moments_retrieve_at_time(IxnParamTraj::ASLEEP,_t_opt,fnames_to_use[_t_opt].size());
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
				itbf->update(_n_t_soln, _time.delta(), dopt, options.local_decay, options.local_decay_factor);
			};
		};

		// Write the basis funcs one last time
		if (options.write) {
			write_bfs(options.dir_write+"F/",n_opt+options.opt_idx_start_writing);
		};
	};



	/********************
	Solve over varying initial conditions
	********************/

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

		// Optimization step translated by the offset
		int i_opt_translated;

		// Iterate over optimization steps
		for (int i_opt=0; i_opt<n_opt; i_opt++)
		{
			std::cout << "Opt step " << i_opt << " / " << n_opt-1 << std::endl;

			// Offset 
			i_opt_translated = i_opt + options.opt_idx_start_writing;

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

				solve_ixn_param_traj();

				// Write
				if (options.write && fnames_to_use[i_opt][i_batch].write) {
					write_ixn_params(options.dir_write+"ixn_params/",i_opt_translated,fnames_to_use[i_opt][i_batch].idxs);
				};

				if (DIAG_SOLVE) { std::cout << "OK" << std::endl; };

				/*****
				Step 2.3 - Solve the variational problem traj
				*****/

				if (DIAG_SOLVE) { std::cout << "Solving var term" << std::endl; };

				solve_var_traj();

				// Write
				if (options.write_var_terms) {
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

				if (options.write && fnames_to_use[i_opt][i_batch].write) {
					write_moments(options.dir_write+"moments/",i_opt_translated,fnames_to_use[i_opt][i_batch].idxs);
				};

				/*****
				Step 2.7 - Gather the update (but dont commit)
				*****/

				for (auto itbf=_bfs.begin(); itbf!=_bfs.end(); itbf++) {
					itbf->update_gather(_n_t_soln, _time.delta(), dopt, options.local_decay, options.local_decay_factor);
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
	OptProblem::OptProblem(std::vector<Dim> dims, std::vector<std::string> species_visible, std::vector<std::string> species_hidden, double t_max, int n_t, int box_length, int lattice_dim) : _impl(new Impl(dims,species_visible,species_hidden,t_max,n_t,box_length,lattice_dim)) {};
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

	void OptProblem::solve_ixn_param_traj() {
		_impl->solve_ixn_param_traj();
	};
	void OptProblem::solve_var_traj() {
		_impl->solve_var_traj();
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





