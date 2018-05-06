#include "../include/bmla.hpp"
#include <iostream>
#include <fstream>
#include "math.h"
#include <sstream>
#include <iomanip>
#include "ixn_param.hpp"
#include "../include/general.hpp"
#include <set>

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

		// Guess
		double _guess;

		// Constructor helpers
		void _clean_up();
		void _copy(const Impl& other);
		void _reset();
		void _shared_constructor(std::string name, DimType type, double guess);

	public:

		// Constructor
		Impl(std::string name, DimType type, double guess);
		Impl(std::string name, DimType type, std::string species, double guess);
		Impl(std::string name, DimType type, std::vector<std::string> species, double guess);
		Impl(std::string name, DimType type, std::vector<std::vector<std::string>> species, double guess);
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

		// Guess
		double guess() const;

		// Get species
		std::vector<std::string> get_species_h() const;
		std::vector<std::string> get_species_b() const;
		std::vector<std::vector<std::string>> get_species_J() const;
		std::vector<std::vector<std::string>> get_species_K() const;
		std::vector<std::vector<std::string>> get_species_W() const;
		std::vector<std::vector<std::string>> get_species_Q() const;

		/********************
		Setters
		********************/

		// Add species
		void add_species_h(std::string species);
		void add_species_b(std::string species);
		void add_species_J(std::string species1, std::string species2);
		void add_species_K(std::string species1, std::string species2, std::string species3);
		void add_species_W(std::string species_visible, std::string species_hidden);
		void add_species_Q(std::string species1, std::string species2, std::string species3, std::string species4);
	};

	/********************
	Constructor
	********************/

	Dim::Impl::Impl(std::string name, DimType type, double guess) {
		_shared_constructor(name,type,guess);
		// Yes any species
		_any_species = true;
	};
	Dim::Impl::Impl(std::string name, DimType type, std::string species, double guess) {
		_shared_constructor(name,type,guess);
		// Not any species
		_any_species = false;
		// Add
		if (_type == H) {
			add_species_h(species);
		} else if (_type == B) {
			add_species_b(species);
		};
	};
	Dim::Impl::Impl(std::string name, DimType type, std::vector<std::string> species, double guess) {
		_shared_constructor(name,type,guess);
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
		} else if (type == K) {
			if (species.size() != 3) {
				std::cerr << "Error! must be 3 species for K" << std::endl;
				exit(EXIT_FAILURE);
			};
			add_species_K(species[0],species[1],species[2]);
		} else if (type == W) {
			if (species.size() != 2) {
				std::cerr << "Error! must be 2 species for W" << std::endl;
				exit(EXIT_FAILURE);
			};
			add_species_W(species[0],species[1]);
 		} else if (type == Q) {
			if (species.size() != 4) {
				std::cerr << "Error! must be 4 species for Q" << std::endl;
				exit(EXIT_FAILURE);
			};
			add_species_Q(species[0],species[1],species[2],species[3]);
		};
	};
	Dim::Impl::Impl(std::string name, DimType type, std::vector<std::vector<std::string>> species, double guess) {
		// Check type
		if (type == B || type != H) {
			std::cerr << "Error! Only for J or K or W or Q." << std::endl;
			exit(EXIT_FAILURE);
		};
		_shared_constructor(name,type,guess);
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
		} else if (type == K) {
			for (auto s_triplet: species) {
				if (s_triplet.size() != 3) {
					std::cerr << "Error! must be 3 species for K" << std::endl;
					exit(EXIT_FAILURE);
				};
				add_species_K(s_triplet[0],s_triplet[1],s_triplet[2]);
			};
		} else if (type == W) {
			for (auto s_pair: species) {
				if (s_pair.size() != 2) {
					std::cerr << "Error! must be 2 species for W" << std::endl;
					exit(EXIT_FAILURE);
				};
				add_species_W(s_pair[0],s_pair[1]);
			};
		} else if (type == Q) {
			for (auto s_quartic: species) {
				if (s_quartic.size() != 4) {
					std::cerr << "Error! must be 4 species for Q" << std::endl;
					exit(EXIT_FAILURE);
				};
				add_species_Q(s_quartic[0],s_quartic[1],s_quartic[2],s_quartic[3]);
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
		_name = other._name;
		_type = other._type;
		_any_species = other._any_species;
		_species = other._species;
		_species_multiple = other._species_multiple;
		_guess = other._guess;
	};
	void Dim::Impl::_reset() {
		_name = "";
		_any_species = false;
		_species.clear();
		_species_multiple.clear();
		_guess = 0.;
	};
	void Dim::Impl::_shared_constructor(std::string name, DimType type, double guess) {
		_name = name;
		_type = type;
		_guess = guess;
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

	// Guess
	double Dim::Impl::guess() const {
		return _guess;
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
	std::vector<std::vector<std::string>> Dim::Impl::get_species_K() const {
		if (_type != K) {
			std::cerr << "Error! Requested species but not of type K." << std::endl;
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
	std::vector<std::vector<std::string>> Dim::Impl::get_species_Q() const {
		if (_type != Q) {
			std::cerr << "Error! Requested species but not of type Q." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _species_multiple;
	};

	/********************
	Setters
	********************/

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
	void Dim::Impl::add_species_K(std::string species1, std::string species2, std::string species3) {
		// Check type
		if (_type != K) {
			std::cerr << "Error! Not K." << std::endl;
			exit(EXIT_FAILURE);
		};
		_any_species = false;
		_species_multiple.push_back(std::vector<std::string>({species1,species2,species3}));
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
	void Dim::Impl::add_species_Q(std::string species1, std::string species2, std::string species3, std::string species4) {
		// Check type
		if (_type != Q) {
			std::cerr << "Error! Not Q." << std::endl;
			exit(EXIT_FAILURE);
		};
		_any_species = false;
		_species_multiple.push_back(std::vector<std::string>({species1,species2,species3,species4}));
	};





















































	/****************************************
	Dim - Impl forwards
	****************************************/

	Dim::Dim(std::string name, DimType type, double guess) : _impl(new Impl(name,type,guess)) {};
	Dim::Dim(std::string name, DimType type, std::string species, double guess) : _impl(new Impl(name,type,species,guess)) {};
	Dim::Dim(std::string name, DimType type, std::vector<std::string> species, double guess) : _impl(new Impl(name,type,species,guess)) {};
	Dim::Dim(std::string name, DimType type, std::vector<std::vector<std::string>> species, double guess) : _impl(new Impl(name,type,species,guess)) {};
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

	// Guess
	double Dim::guess() const {
		return _impl->guess();
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
	std::vector<std::vector<std::string>> Dim::get_species_K() const {
		return _impl->get_species_K();
	};
	std::vector<std::vector<std::string>> Dim::get_species_W() const {
		return _impl->get_species_W();
	};
	std::vector<std::vector<std::string>> Dim::get_species_Q() const {
		return _impl->get_species_Q();
	};

	/********************
	Setters
	********************/

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
	void Dim::add_species_K(std::string species1, std::string species2, std::string species3) {
		_impl->add_species_K(species1,species2,species3);
	};
	void Dim::add_species_W(std::string species_visible, std::string species_hidden) {
		_impl->add_species_W(species_visible,species_hidden);
	};
	void Dim::add_species_Q(std::string species1, std::string species2, std::string species3, std::string species4) {
		_impl->add_species_Q(species1,species2,species3,species4);
	};










































	/****************************************
	BMLA - IMPLEMENTATION
	****************************************/

	class BMLA::Impl {

	private:

		// Number of dimensions
		int _n_param;

		// List of interaction parameters
		std::list<IxnParam> _ixn_params;

		// Species present
		std::list<Species> _species;
		std::list<HiddenSpecies> _hidden_species;

		// List of hidden units, and flag if they exist
		bool _hidden_layer_exists;
		std::list<HiddenUnit> _hidden_units;
		std::list<ConnectionVH> _conn_vh;

		// Counters
		std::list<Counter> _counters;

		// The lattice to learn
		Lattice _latt;

		/********************
		Private methods
		********************/

		// Print
		void _print_ixn_params(bool new_line=true) const;
		void _print_moments(bool new_line=true) const;
		void _print_mse(bool new_line=true) const;

		// Get the mse
		double _get_mse() const;

		// Add a hidden unit		
		void _add_hidden_unit(std::vector<std::string> species_possible, std::vector<Site*> conn_sites, std::vector<std::string> w_params, std::vector<std::string> b_params);
		std::vector<Site*> _get_sites(std::vector<int> &lattice_idxs);
		std::vector<Site*> _get_sites(std::vector<std::vector<int>> &lattice_idxs);

		// Search functions
		Species* _not_nullptr(Species* ptr);
		HiddenSpecies* _not_nullptr(HiddenSpecies *ptr);
		IxnParam* _not_nullptr(IxnParam* ptr);
		Counter* _not_nullptr(Counter* ptr);
		// Find species
		Species* _find_species(std::string name);
		// Find hidden species
		HiddenSpecies* _find_hidden_species(std::string name);
		// Find Ixn param
		IxnParam* _find_ixn_param_by_name(std::string name);
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
		void _copy(const Impl& other);
		void _reset();

	public:

		// Constructor
		Impl(std::vector<Dim> dims, std::vector<std::string> species_visible, std::vector<std::string> species_hidden, int box_length, int lattice_dim);
		Impl(const Impl& other);
		Impl(Impl&& other);
		Impl& operator=(const Impl& other);
	    Impl& operator=(Impl&& other);
		~Impl();

		// Set a parameter for dim
		void set_param_for_dim(std::string dim_name, double val);

		// Set a fixed moment
		void set_fixed_awake_moment_for_dim(std::string dim_name, double val);

		// Any dim
		void add_hidden_unit(std::vector<std::string> species_possible, std::vector<std::vector<int>> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params);
		// 1D specific
		void add_hidden_unit(std::vector<std::string> species_possible, std::vector<int> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params);
		// Validate
		void validate_hidden() const;

		// Solve for the h,j corresponding to a given lattice
		void solve(std::vector<std::string> fnames, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolveBMLA options);

		// At the current ixns params, sample and report the specified moments
		void sample(int batch_size, int n_cd_steps, OptionsSampling options);
		// Internal
		void _sample(int n_cd_steps, OptionsSampling options);

		// Update the initial params
		void read(std::string fname);

		// Write out the solutions
		void write(std::string fname, bool write_idx_opt_step, int idx_opt_step, bool append);
		void write_ave(std::string fname, int last_n_steps, bool write_idx_opt_step, int idx_opt_step, bool append);
		void write_moments(std::string fname, int idx_opt_step, bool append);
	};

































	/****************************************
	BMLA - IMPLEMENTATION DEFINITIONS
	****************************************/

	/********************
	Constructor
	********************/

	BMLA::Impl::Impl(std::vector<Dim> dims, std::vector<std::string> species_visible, std::vector<std::string> species_hidden, int box_length, int lattice_dim) : _latt(lattice_dim,box_length) {

		// Number of dims
		_n_param = dims.size();
		// No hidden layers yet
		_hidden_layer_exists = false;

		// Fix dims that have any species specification
		for (auto d=dims.begin(); d!=dims.end(); d++) {
			if (d->any_species()) {
				if (d->type()==H) {
					for (auto s: species_visible) {
						d->add_species_h(s);
					};
				} else if (d->type()==B) {
					for (auto s: species_hidden) {
						d->add_species_b(s);
					};
				} else if (d->type()==J) {
					for (auto s1: species_visible) {
						for (auto s2: species_visible) {
							d->add_species_J(s1,s2);
						};
					};
				} else if (d->type()==K) {
					for (auto s1: species_visible) {
						for (auto s2: species_visible) {
							for (auto s3: species_visible) {
								d->add_species_K(s1,s2,s3);
							};
						};
					};
				} else if (d->type()==W) {
					for (auto sv: species_visible) {
						for (auto sh: species_hidden) {
							d->add_species_W(sv,sh);
						};
					};
				} else if (d->type()==Q) {
					for (auto s1: species_visible) {
						for (auto s2: species_visible) {
							for (auto s3: species_visible) {
								for (auto s4: species_visible) {
									d->add_species_Q(s1,s2,s3,s4);
								};
							};
						};
					};
				};
			};
		};

		// Check if hidden layer exists
		for (auto d=dims.begin(); d!=dims.end(); d++) {
			if (d->type()==B || d->type()==W) {
				_hidden_layer_exists = true;
				break;
			};
		};

		// Create the visible species and add to the lattice
		for (auto s: species_visible) {
			// Make species
			_species.push_back(Species(s));

			// Add to the lattice
			_latt.add_species_possibility(&(_species.back()));
		};

		// Create the hidden species
		for (auto s: species_hidden) {
			// Make species
			_hidden_species.push_back(HiddenSpecies(s));
		};

		// Create the interaction params
		Species *sp,*sp1,*sp2,*sp3,*sp4;
		HiddenSpecies *sph;
		for (auto d=dims.begin(); d!=dims.end(); d++) {

			if (d->type()==H) {

				// Create
				_ixn_params.push_back(IxnParam(d->name(),Hp,d->guess()));

				for (auto s: d->get_species_h()) {
					// Find and add the species
					sp = _not_nullptr(_find_species(s));
					_ixn_params.back().add_species(sp);
					// Add ixn to the species
					sp->add_h_ptr(&_ixn_params.back());
					// Make sure we have a counter
					_add_counter(s);	
				};

			} else if (d->type()==J) { 

				// Create
				_ixn_params.push_back(IxnParam(d->name(),Jp,d->guess()));

				for (auto s: d->get_species_J()) {
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

			} else if (d->type()==K) {

				// Check that Lattice dim is 1 - only 1 is allowed!
				if (lattice_dim != 1) {
					std::cerr << "Error: triplets are currently only supported for lattice of dim 1" << std::endl;
					exit(EXIT_FAILURE);
				};

				// Create
				_ixn_params.push_back(IxnParam(d->name(),Kp,d->guess()));

				for (auto s: d->get_species_K()) {
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

			} else if (d->type()==Q) {

				// Check that Lattice dim is 1 - only 1 is allowed!
				if (lattice_dim != 1) {
					std::cerr << "Error: quartics are currently only supported for lattice of dim 1" << std::endl;
					exit(EXIT_FAILURE);
				};

				// Create
				_ixn_params.push_back(IxnParam(d->name(),Qp,d->guess()));

				for (auto s: d->get_species_Q()) {
					// Find and add the species
					sp1 = _not_nullptr(_find_species(s[0]));
					sp2 = _not_nullptr(_find_species(s[1]));
					sp3 = _not_nullptr(_find_species(s[2]));
					sp4 = _not_nullptr(_find_species(s[3]));
					_ixn_params.back().add_species(sp1,sp2,sp3,sp4);
					// Add ixn to the species
					sp1->add_q_ptr(sp2,sp3,sp4,&_ixn_params.back());
					sp2->add_q_ptr(sp1,sp3,sp4,&_ixn_params.back());
					sp3->add_q_ptr(sp1,sp2,sp4,&_ixn_params.back());
					sp4->add_q_ptr(sp1,sp2,sp3,&_ixn_params.back());
					// Make sure we have a counter
					_add_counter(s[0],s[1],s[2],s[3]);			
				};

				// Init structure for lattice
				_latt.init_quartic_structure();

			} else if (d->type()==W) { 

				// Create
				_ixn_params.push_back(IxnParam(d->name(),Wp,d->guess()));

				for (auto s: d->get_species_W()) {
					// Find and add the species
					sp = _not_nullptr(_find_species(s[0]));
					sph = _not_nullptr(_find_hidden_species(s[1]));
					_ixn_params.back().add_species(sp,sph);

					// No need to add to species - handled by ConnectionVH class
				};

			} else if (d->type()==B) {

				// Create
				_ixn_params.push_back(IxnParam(d->name(),Bp,d->guess()));

				for (auto s: d->get_species_b()) {
					// Find and add the species
					sph = _not_nullptr(_find_hidden_species(s));
					_ixn_params.back().add_species(sph);

					// No need to add to species - handled by ConnectionVH class
				};

			};
		};
	};
	BMLA::Impl::Impl(const Impl& other) : _latt(other._latt) {
		_copy(other);
	};
	BMLA::Impl::Impl(Impl&& other) : _latt(other._latt) {
		_copy(other);
		other._reset();
	};
    BMLA::Impl& BMLA::Impl::operator=(const Impl& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    BMLA::Impl& BMLA::Impl::operator=(Impl&& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
    };
	BMLA::Impl::~Impl()
	{
		_clean_up();
	};

	void BMLA::Impl::_clean_up() {
		// Nothing...
	};
	void BMLA::Impl::_copy(const Impl& other) {
		_n_param = other._n_param;
		_ixn_params = other._ixn_params;
		_species = other._species;
		_hidden_layer_exists = other._hidden_layer_exists;
		_hidden_units = other._hidden_units;
		_conn_vh = other._conn_vh;
		_counters = other._counters;
		_latt = other._latt;
	};
	void BMLA::Impl::_reset() {
		_n_param = 0;
		_ixn_params.clear();
		_species.clear();
		_hidden_layer_exists = false;
		_hidden_units.clear();
		_conn_vh.clear();
		_counters.clear();
		_latt = Lattice(0,0);
	};

	/********************
	Print (compare) moments
	********************/

	void BMLA::Impl::_print_ixn_params(bool new_line) const {
		std::cout << "   Ixn params: ";
 		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			std::cout << it->get() << " ";
		};
		if (new_line) {
			std::cout << std::endl;
		};
	};

	void BMLA::Impl::_print_moments(bool new_line) const {
		std::cout << "   Moments: ";
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			std::cout << it->get_moment(IxnParam::ASLEEP) << " (" << it->get_moment(IxnParam::AWAKE) << "), ";
		};
		if (new_line) {
			std::cout << std::endl;
		};
	};
	void BMLA::Impl::_print_mse(bool new_line) const {
		std::cout << "   MSE: " << _get_mse() << "%";
		if (new_line) {
			std::cout << std::endl;
		};
	};

	/********************
	Get MSE
	********************/

	double BMLA::Impl::_get_mse() const {
		double mse=0.0;
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			mse += abs(it->moments_diff())/abs(it->get_moment(IxnParam::AWAKE));
		};
		return 100*mse/_n_param;
	};

	/********************
	Add counter for some species or nns
	********************/

	void BMLA::Impl::_add_counter(std::string s) {
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
	void BMLA::Impl::_add_counter(std::string s1, std::string s2) {
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
	void BMLA::Impl::_add_counter(std::string s1, std::string s2, std::string s3) {
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
	void BMLA::Impl::_add_counter(std::string s1, std::string s2, std::string s3, std::string s4) {
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

	/********************
	Set a parameter for dim
	********************/

	void BMLA::Impl::set_param_for_dim(std::string dim_name, double val) {
		// Find
		IxnParam *ip = _not_nullptr(_find_ixn_param_by_name(dim_name));

		// Guaranteed not to be null
		ip->set_val(val);
	};

	/********************
	Set a fixed value for a moment to learn
	********************/

	void BMLA::Impl::set_fixed_awake_moment_for_dim(std::string dim_name, double val) {
		// Find
		IxnParam *ip = _not_nullptr(_find_ixn_param_by_name(dim_name));

		// Guaranteed not to be null
		ip->set_fixed_awake_moment(val);
	};

	/********************
	Find hidden unit connections
	********************/

	// Any dim
	std::vector<Site*> BMLA::Impl::_get_sites(std::vector<std::vector<int>> &lattice_idxs) {
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
	std::vector<Site*> BMLA::Impl::_get_sites(std::vector<int> &lattice_idxs) {
		std::vector<Site*> conns;
		for (auto c: lattice_idxs) {
			conns.push_back(_latt.get_site(c));
		};
		return conns;
	};

	/********************
	Add hidden unit
	********************/

	void BMLA::Impl::add_hidden_unit(std::vector<std::string> species_possible, std::vector<std::vector<int>> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params) {
		// Find sites indicated by connections
		std::vector<Site*> conns = _get_sites(lattice_idxs);
		// Make
		_add_hidden_unit(species_possible, conns, w_params, b_params);
	};
	void BMLA::Impl::add_hidden_unit(std::vector<std::string> species_possible, std::vector<int> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params) {
		// Find sites indicated by connections
		std::vector<Site*> conns = _get_sites(lattice_idxs);
		// Make
		_add_hidden_unit(species_possible, conns, w_params, b_params);
	};

	/********************
	Add hidden unit internal
	********************/

	void BMLA::Impl::_add_hidden_unit(std::vector<std::string> species_possible, std::vector<Site*> conn_sites, std::vector<std::string> w_params, std::vector<std::string> b_params) {

		// Find the ixn params listed
		std::vector<IxnParam*> ip_w;
		std::vector<IxnParam*> ip_b;
		IxnParam *ip;
		for (auto w: w_params) {
			ip = _not_nullptr(_find_ixn_param_by_name(w));
			ip_w.push_back(ip);
		};
		for (auto b: b_params) {
			ip = _not_nullptr(_find_ixn_param_by_name(b));
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
	Validate hidden layer
	********************/

	void BMLA::Impl::validate_hidden() const {
		_latt.validate_graph();
	};

	/********************
	Solve for the h,j corresponding to a given lattice
	********************/

	void BMLA::Impl::solve(std::vector<std::string> fnames, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolveBMLA options)
	{
		// Check size of filenames
		if (options.use_single_lattice == true && fnames.size() != 1) {
			std::cerr << "Error! In _use_single_lattice mode, only provide one filename!" << std::endl;
			exit(EXIT_FAILURE); 
		};

		// Reset the params to the guesses, and the moments to 0
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			it->reset_to_guess();
			it->moments_reset(IxnParam::AWAKE);
			it->moments_reset(IxnParam::ASLEEP);

			// Reset solution trajs as needed
			if (options.track_soln_traj) {
				it->set_track_soln_traj(true);
				it->reset_soln_traj();
			};
		};

		// Opt step with offset
		int i_opt = options.opt_idx_start_writing;

		// Write the initial point for the solution
		if (options.write_soln_traj) {
			write(options.fname_write_soln_traj,true,i_opt,options.append);
		};

		// Iterate over optimization steps
		for (int i_opt_from_zero=0; i_opt_from_zero<n_opt; i_opt_from_zero++)
		{
			if (options.verbose) {
				std::cout << "Opt step: " << i_opt_from_zero << " / " << n_opt << std::endl;
			};

			// Offset
			i_opt = i_opt_from_zero + options.opt_idx_start_writing;

			// Check MSE to see if quit
			if (options.mse_quit_mode) {
				if (_get_mse() < options.mse_quit) {
					// Quit!
					std::cout << "--- MSE is low enough: " << _get_mse() << " < " << options.mse_quit << " quitting! ---" << std::endl;
					break;
				};
			};

			if (options.nesterov) {
				// If first opt step, do nothing, but set the "prev" point to the current to initialize
				if (i_opt_from_zero == 0) {
					for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
						it->nesterov_set_prev_equal_curr();
					};
				} else {
					// Move to the intermediate point
					for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
						it->nesterov_move_to_intermediate_pt(i_opt);
					};
				};
			};

			// Reset all moments
			for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
				it->moments_reset(IxnParam::AWAKE);
				it->moments_reset(IxnParam::ASLEEP);
			};

			// Go over the batch
			if (options.verbose) { std::cout << "   " << std::flush; };
			for (int i_batch=0; i_batch<batch_size; i_batch++)
			{
				if (options.verbose) {
					std::cout << "." << std::flush;
				};

				// Reset the lattice by reading in a random
				if (!options.start_CD_with_empty_latt)
				{	
					_latt.read_from_file(fnames[randI(0,fnames.size()-1)],options.awake_visible_are_binary);
				} else {
					// Don't read in, just start from empty
					_latt.clear();
				};

				// Activate hidden
				if (_hidden_layer_exists) {
					for (auto ithu = _hidden_units.begin(); ithu != _hidden_units.end(); ithu++) {
						if (options.awake_hidden_are_binary) {
							// Binary
							ithu->activate(true);
						} else {
							// Probabilistic
							ithu->activate(false);
						};
					};
				};

				// Record awake moments
				for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
					it->moments_retrieve(IxnParam::AWAKE, batch_size);
				};

				// Do the contrastive divergence
				for (int cd_step=0; cd_step<n_cd_steps; cd_step++)
				{
					// Sample
					if (cd_step != n_cd_steps-1) {
						_latt.sample(options.asleep_visible_are_binary);
					} else {
						_latt.sample(options.asleep_final_visible_are_binary);
					};

					// Activate hidden
					if (_hidden_layer_exists) {
						for (auto ithu = _hidden_units.begin(); ithu != _hidden_units.end(); ithu++) {
							if (cd_step != n_cd_steps-1) {
								ithu->activate(options.asleep_hidden_are_binary);
							} else {
								ithu->activate(options.asleep_final_hidden_are_binary);
							};
						};
					};

					// Hidden layer:
					/*
					std::cout << "hidden layer:" << std::endl;
					for (auto ithu = _hidden_units.begin(); ithu != _hidden_units.end(); ithu++) {
						ithu->print_conns(false);
						std::cout << ithu->get() << std::endl;
					};
					*/
				};

				// Record asleep moments
				//std::cout << "---- retrieving asleep moments ----" << std::endl;
				for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
					it->moments_retrieve(IxnParam::ASLEEP, batch_size);
				};
				//std::cout << "------------------------------" << std::endl;

			};

			// Print out the MSE
			if (options.verbose) {
				_print_ixn_params(false);
				_print_moments(false);
				_print_mse();
			};

			// Write the moments if needed
			if (options.write_moment_traj) {
				if (i_opt_from_zero==0) {
					// Possibly make new depending on append flag
					write_moments(options.fname_write_moment_traj,i_opt+1,options.append);
				} else {
					// Append
					write_moments(options.fname_write_moment_traj,i_opt+1,true);
				};
			};

			// Update the params
			for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
				it->update(dopt,options.l2_reg_mode,options.l2_lambda);
			};

			// Write the new solution (append)
			if (options.write_soln_traj) {
				write(options.fname_write_soln_traj,true,i_opt+1,true);
			};
		};

		// Report final
		std::cout << "--- Final ---" << std::endl;
		_print_ixn_params();
		_print_moments();
		_print_mse();
	};

	/********************
	At the current ixns params, sample and report the specified moments
	********************/

	// Internal sampling function
	void BMLA::Impl::_sample(int n_cd_steps, OptionsSampling options) {

		// Reset stored counters
		for (auto it_ctr = _counters.begin(); it_ctr != _counters.end(); it_ctr++) {
			if (it_ctr->report_during_sampling()) {
				it_ctr->storage_clear();
			};
		};

		// Start by populating lattice randomly
		_latt.populate_randomly();

		// Activate hidden
		if (_hidden_layer_exists) {
			for (auto ithu = _hidden_units.begin(); ithu != _hidden_units.end(); ithu++) {
				if (options.awake_hidden_are_binary) {
					ithu->activate(true);
				} else {
					ithu->activate(false);
				};
			};
		};

		// Store counts
		for (auto it_ctr = _counters.begin(); it_ctr != _counters.end(); it_ctr++) {
			if (it_ctr->report_during_sampling()) {
				it_ctr->storage_committ_current_count();
			};
		};

		// CD
		for (int cd_step=0; cd_step<n_cd_steps; cd_step++) {

			// Sample visibles
			if (cd_step != n_cd_steps-1) {
				_latt.sample(options.asleep_visible_are_binary);
			} else {
				_latt.sample(options.asleep_final_visible_are_binary);
			};

			// Activate hidden
			if (_hidden_layer_exists) {
				for (auto ithu = _hidden_units.begin(); ithu != _hidden_units.end(); ithu++) {
					if (cd_step != n_cd_steps-1) {
						ithu->activate(options.asleep_hidden_are_binary);
					} else {
						ithu->activate(options.asleep_final_hidden_are_binary);
					}
				};
			};

			// Store counts
			for (auto it_ctr = _counters.begin(); it_ctr != _counters.end(); it_ctr++) {
				if (it_ctr->report_during_sampling()) {
					it_ctr->storage_committ_current_count();
				};
			};
		};
	};

	void BMLA::Impl::sample(int batch_size, int n_cd_steps, OptionsSampling options) {

		// Check writing
		if (options.write_traj && options.fname_write_traj == "") {
			std::cerr << "ERROR: provide a filename for writing!" << std::endl;
			exit(EXIT_FAILURE);
		};

		// Create any counters, as needed
		Counter *ctr;
		for (auto sp: options.report_counts) {
			_add_counter(sp);
			ctr = _not_nullptr(_find_ctr_by_species(sp));
			ctr->set_report_during_sampling(true);
		};
		for (auto sp_nbrs: options.report_nns) {
			_latt.init_nn_structure();
			_add_counter(sp_nbrs.first,sp_nbrs.second);
			ctr = _not_nullptr(_find_ctr_by_species(sp_nbrs.first,sp_nbrs.second));
			ctr->set_report_during_sampling(true);
		};
		for (auto sp_triplets: options.report_triplets) {
			_latt.init_triplet_structure();
			_add_counter(sp_triplets[0],sp_triplets[1],sp_triplets[2]);
			ctr = _not_nullptr(_find_ctr_by_species(sp_triplets[0],sp_triplets[1],sp_triplets[2]));
			ctr->set_report_during_sampling(true);
		};
		for (auto sp_quartics: options.report_quartics) {
			_latt.init_quartic_structure();
			_add_counter(sp_quartics[0],sp_quartics[1],sp_quartics[2],sp_quartics[3]);
			ctr = _not_nullptr(_find_ctr_by_species(sp_quartics[0],sp_quartics[1],sp_quartics[2],sp_quartics[3]));
			ctr->set_report_during_sampling(true);
		};

		// Reset stored counters
		for (auto it_ctr = _counters.begin(); it_ctr != _counters.end(); it_ctr++) {
			if (it_ctr->report_during_sampling()) {
				it_ctr->storage_averaged_clear();
			};
		};		

		// Go through batches
		for (int batch_no=0; batch_no<batch_size; batch_no++)
		{
			std::cout << "." << std::flush;

			// Sample
			_sample(n_cd_steps,options);

			// Store the moments
			for (auto it_ctr = _counters.begin(); it_ctr != _counters.end(); it_ctr++) {
				if (it_ctr->report_during_sampling()) {
					it_ctr->storage_averaged_committ_current_traj(batch_size);
				};
			};

		};
		std::cout << std::endl;

		// Report the moments
		if (options.verbose)
		{
			for (auto it_ctr = _counters.begin(); it_ctr != _counters.end(); it_ctr++) {
				if (it_ctr->report_during_sampling()) {
					it_ctr->storage_averaged_print();
				};
			};
		};

		// Write if needed
		if (options.write_traj) {
			std::ofstream f;
			f.open(options.fname_write_traj);
			for (auto it_ctr = _counters.begin(); it_ctr != _counters.end(); it_ctr++) {
				if (it_ctr->report_during_sampling()) {
					it_ctr->storage_averaged_write(f);
				};
			};
			f.close();
		};
	};

	/********************
	Read initial guess
	********************/

	void BMLA::Impl::read(std::string fname) 
	{
		std::ifstream f;
		f.open(fname);
		std::string ixn_name="";
		std::string guess="";
		std::string line;
		std::istringstream iss;
		IxnParam* ip;
		if (f.is_open()) { // make sure we found it
			while (getline(f,line)) {
				if (line == "") { continue; };
				iss = std::istringstream(line);
			    iss >> ixn_name;
			    iss >> guess;
		    	// Add
		    	ip = _not_nullptr(_find_ixn_param_by_name(ixn_name));
			    ip->set_guess(atof(guess.c_str()));
			    ip->set_val(atof(guess.c_str()));
		    	ixn_name=""; guess="";
			};
		};
		f.close();
	};

	/********************
	Write the solutions
	********************/

	void BMLA::Impl::write(std::string fname, bool write_idx_opt_step, int idx_opt_step, bool append) {
		std::ofstream f;
		if (append) {
			f.open(fname, std::ofstream::out | std::ofstream::app);
		} else {
			f.open(fname);
		};
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			if (write_idx_opt_step) {
				f << idx_opt_step << " " << it->name() << " " << it->get() << "\n";
			} else {
				f << it->name() << " " << it->get() << "\n";
			};
		};
		f.close();	
	};
	void BMLA::Impl::write_ave(std::string fname, int last_n_steps, bool write_idx_opt_step, int idx_opt_step, bool append) {
		std::ofstream f;
		if (append) {
			f.open(fname, std::ofstream::out | std::ofstream::app);
		} else {
			f.open(fname);
		};

		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			if (write_idx_opt_step) {
				f << idx_opt_step << " " << it->name() << " " << it->get_ave(last_n_steps) << "\n";
			} else {
				f << it->name() << " " << it->get() << "\n";
			};
		};
		f.close();	
	};
	void BMLA::Impl::write_moments(std::string fname, int idx_opt_step, bool append) {
		std::ofstream f;
		if (append) {
			f.open(fname, std::ofstream::out | std::ofstream::app);
		} else {
			f.open(fname);
		};
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			f << idx_opt_step << " " << it->name() << " " << it->get_moment(IxnParam::AWAKE) << " " << it->get_moment(IxnParam::ASLEEP) << "\n";
		};
		f.close();	
	};

	/********************
	Search functions
	********************/

	Species* BMLA::Impl::_not_nullptr(Species* ptr) {
		if (!ptr) {
			std::cerr << "ERROR: could not find species!" << std::endl;
			exit(EXIT_FAILURE);
		};
		return ptr;
	};
	HiddenSpecies* BMLA::Impl::_not_nullptr(HiddenSpecies* ptr) {
		if (!ptr) {
			std::cerr << "ERROR: could not find hidden species!" << std::endl;
			exit(EXIT_FAILURE);
		};
		return ptr;
	};
	IxnParam* BMLA::Impl::_not_nullptr(IxnParam* ptr) {
		if (!ptr) {
			std::cerr << "ERROR: could not find IxnParam!" << std::endl;
			exit(EXIT_FAILURE);
		};
		return ptr;
	};
	Counter* BMLA::Impl::_not_nullptr(Counter* ptr) {
		if (!ptr) {
			std::cerr << "ERROR: could not find counter!" << std::endl;
			exit(EXIT_FAILURE);
		};
		return ptr;
	};
	Species* BMLA::Impl::_find_species(std::string name) {
		for (auto it=_species.begin(); it!=_species.end(); it++) {
			if (it->name() == name) {
				return &*it;
			};
		};
		return nullptr;
	};
	HiddenSpecies* BMLA::Impl::_find_hidden_species(std::string name) {
		for (auto it=_hidden_species.begin(); it!=_hidden_species.end(); it++) {
			if (it->name() == name) {
				return &*it;
			};
		};
		return nullptr;
	};
	IxnParam* BMLA::Impl::_find_ixn_param_by_name(std::string name) {
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			if (it->name() == name) {
				return &*it;
			};
		};
		return nullptr;
	};
	Counter* BMLA::Impl::_find_ctr_by_species(std::string s) {
		for (auto it=_counters.begin(); it!=_counters.end(); it++) {
			if (it->is_counting_species(s)) {
				return &*it;
			};
		};
		return nullptr;
	};
	Counter* BMLA::Impl::_find_ctr_by_species(std::string s1, std::string s2) {
		for (auto it=_counters.begin(); it!=_counters.end(); it++) {
			if (it->is_counting_species(s1,s2)) {
				return &*it;
			};
		};
		return nullptr;
	};
	Counter* BMLA::Impl::_find_ctr_by_species(std::string s1, std::string s2, std::string s3) {
		for (auto it=_counters.begin(); it!=_counters.end(); it++) {
			if (it->is_counting_species(s1,s2,s3)) {
				return &*it;
			};
		};
		return nullptr;
	};
	Counter* BMLA::Impl::_find_ctr_by_species(std::string s1, std::string s2, std::string s3, std::string s4) {
		for (auto it=_counters.begin(); it!=_counters.end(); it++) {
			if (it->is_counting_species(s1,s2,s3,s4)) {
				return &*it;
			};
		};
		return nullptr;
	};







































	/****************************************
	BMLA IMPL forwards
	****************************************/

	// Constructor
	BMLA::BMLA(std::vector<Dim> dims, std::vector<std::string> species_visible, std::vector<std::string> species_hidden, int box_length, int lattice_dim) : _impl(new Impl(dims,species_visible,species_hidden,box_length,lattice_dim)) {};
	BMLA::BMLA(const BMLA& other) : _impl(new Impl(*other._impl)) {};
	BMLA::BMLA(BMLA&& other) = default;
	BMLA& BMLA::operator=(BMLA other) {
        _impl = std::move(other._impl);
        return *this; 
	};
	BMLA::~BMLA() = default;

	// Set a parameter for dim
	void BMLA::set_param_for_dim(std::string dim_name, double val) {
		_impl->set_param_for_dim(dim_name,val);
	};

	// Set a fixed value for a moment
	void BMLA::set_fixed_awake_moment_for_dim(std::string dim_name, double val) {
		_impl->set_fixed_awake_moment_for_dim(dim_name,val);
	};

	// Add hidden unit for any dim
	void BMLA::add_hidden_unit(std::vector<std::string> species_possible, std::vector<std::vector<int>> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params) {
		_impl->add_hidden_unit(species_possible, lattice_idxs, w_params, b_params);
	};
	// 1D specific
	void BMLA::add_hidden_unit(std::vector<std::string> species_possible, std::vector<int> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params) {
		_impl->add_hidden_unit(species_possible, lattice_idxs,w_params,b_params);
	};
	// Validate
	void BMLA::validate_hidden() const {
		_impl->validate_hidden();
	};

	// Solve for the h,j corresponding to a given lattice
	void BMLA::solve(std::string fname, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolveBMLA options) {
		std::vector<std::string> fnames;
		fnames.push_back(fname);
		_impl->solve(fnames,n_opt,batch_size,n_cd_steps,dopt,options);
	};

	// Solve for the h,j corresponding to a given lattice
	void BMLA::solve(std::vector<std::string> fnames, int n_opt, int batch_size, int n_cd_steps, double dopt, OptionsSolveBMLA options) {
		_impl->solve(fnames,n_opt,batch_size,n_cd_steps,dopt,options);
	};

	// At the current ixns params, sample and report the specified moments
	void BMLA::sample(int batch_size, int n_cd_steps, OptionsSampling options) {
		_impl->sample(batch_size, n_cd_steps, options);
	};

	// Update the initial params
	void BMLA::read(std::string fname) {
		_impl->read(fname);
	};

	// Write out the solutions
	void BMLA::write(std::string fname, bool append) {
		_impl->write(fname,false,0,append);
	};
	void BMLA::write(std::string fname, int idx, bool append) {
		_impl->write(fname,true,idx,append);
	};
	void BMLA::write_ave(std::string fname, int last_n_steps, bool append) {
		_impl->write_ave(fname,last_n_steps,false,0,append);
	};
	void BMLA::write_ave(std::string fname, int last_n_steps, int idx, bool append) {
		_impl->write_ave(fname,last_n_steps,true,idx,append);
	};
};


