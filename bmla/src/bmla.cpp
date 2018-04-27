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
	Struct for specifying a dimension
	****************************************/

	Dim::Dim(std::string name, DimType type, double guess) {
		if (type != B && type != W) {
			std::cerr << "ERROR! Species-less types other than B,W are not supported yet." << std::endl;
			exit(EXIT_FAILURE);
		};
		this->name = name;
		this->type = type;
		this->all_species = true;
		this->species1 = "";
		this->species2 = "";
		this->species3 = "";
		this->guess = guess;
	};
	Dim::Dim(std::string name, DimType type, std::string species, double guess) : Dim(name, type, species, "", "", guess) {};
	Dim::Dim(std::string name, DimType type, std::string species1, std::string species2, double guess) : Dim(name, type, species1, species2, "", guess) {};
	Dim::Dim(std::string name, DimType type, std::string species1, std::string species2, std::string species3, double guess)
	{
		if (type == B || type == H || type == W) {
			if (species1 == "" || species2 != "" || species3 != "") {
				std::cerr << "ERROR! Dim specification is incorrect for H, B or W." << std::endl;
				exit(EXIT_FAILURE);
			};
		} else if (type == J) {
			if (species1 == "" || species2 == "" || species3 != "") {
				std::cerr << "ERROR! Dim specification is incorrect for J." << std::endl;
				exit(EXIT_FAILURE);
			};
		} else if (type == K) {
			if (species1 == "" || species2 == "" || species3 == "") {
				std::cerr << "ERROR! Dim specification is incorrect for K." << std::endl;
				exit(EXIT_FAILURE);
			};
		};

		this->name = name;
		this->type = type;
		this->all_species = false;
		this->species1 = species1;
		this->species2 = species2;
		this->species3 = species3;
		this->guess = guess;
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

		// List of hidden units, and flag if they exist
		bool _hidden_layer_exists;
		std::list<HiddenUnit> _hidden_units;

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
		void _add_hidden_unit(std::vector<Site*> conn_sites, std::vector<std::string> w_params, std::vector<std::string> b_params);
		std::vector<Site*> _get_sites(std::vector<int> &lattice_idxs);
		std::vector<Site*> _get_sites(std::vector<std::vector<int>> &lattice_idxs);

		// Search functions
		Species* _not_nullptr(Species* ptr);
		IxnParam* _not_nullptr(IxnParam* ptr);
		Counter* _not_nullptr(Counter* ptr);
		// Find species
		Species* _find_species(std::string name);
		// Find Ixn param
		IxnParam* _find_ixn_param_by_name(std::string name);
		IxnParam* _find_ixn_param_by_species(IxnParamType type, std::string s);
		IxnParam* _find_ixn_param_by_species(IxnParamType type, std::string s1, std::string s2);	
		IxnParam* _find_ixn_param_by_species(IxnParamType type, std::string s1, std::string s2, std::string s3);
		IxnParam* _find_ixn_param_for_any_species(IxnParamType type);
		// Find counter
		Counter* _find_ctr_by_species(std::string s);
		Counter* _find_ctr_by_species(std::string s1, std::string s2);
		Counter* _find_ctr_by_species(std::string s1, std::string s2, std::string s3);
		Counter* _find_ctr_by_species(std::string s1, std::string s2, std::string s3, std::string s4);

		// Constructor helpers
		void _clean_up();
		void _copy(const Impl& other);
		void _reset();

	public:

		// Constructor
		Impl(std::vector<Dim> dims, int box_length, int lattice_dim);
		Impl(Impl&& other);
	    Impl& operator=(Impl&& other);
		~Impl();

		// Set a parameter for dim
		void set_param_for_dim(std::string dim_name, double val);

		// Any dim
		void add_hidden_unit(std::vector<std::vector<int>> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params);
		// 1D specific
		void add_hidden_unit(std::vector<int> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params);
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

		// Add a counter for some species or nns
		void add_counter(std::string s);
		void add_counter(std::string s1, std::string s2);
		void add_counter(std::string s1, std::string s2, std::string s3);
		void add_counter(std::string s1, std::string s2, std::string s3, std::string s4);
	};

	/****************************************
	BMLA - IMPLEMENTATION DEFINITIONS
	****************************************/

	/********************
	Constructor
	********************/

	BMLA::Impl::Impl(std::vector<Dim> dims, int box_length, int lattice_dim) : _latt(lattice_dim,box_length) {

		// Number of dims
		_n_param = dims.size();
		// No hidden layers yet
		_hidden_layer_exists = false;

		// Tell the lattice about what dims exist for sampling
		for (auto d: dims) {
			if (d.type==H) { 
				_latt.set_sampling_flag_exists_h(true);
			} else if (d.type==J) {
				_latt.set_sampling_flag_exists_j(true);
				// setup lattice for NNs
				_latt.init_nn_structure();
			} else if (d.type==K) {
				_latt.set_sampling_flag_exists_k(true);
				// set up lattice for triplets
				_latt.init_triplet_structure();
			} else if (d.type==W) {
				_latt.set_sampling_flag_exists_w(true);
				// Mark that it exists
				_hidden_layer_exists = true;
			} else if (d.type==B) {
				// No need to tell lattice, since it only affects hidden unit activation and not sampling
				// Mark that hidden layer exists
				_hidden_layer_exists = true;
			};
		};

		// Compile the set of species that exist
		std::set<std::string> species;
		for (auto d: dims) {
			if (d.species1 != "") { species.insert(d.species1); };
			if (d.species2 != "") { species.insert(d.species2); };
			if (d.species3 != "") { species.insert(d.species3); };
		};

		// Create the species and add to the lattice
		std::vector<Species*> sp_all; // all species
		for (auto s: species) {
			// Make species
			_species.push_back(Species(s));

			// Add to the lattice
			_latt.add_species_possibility(&(_species.back()));

			// Keep track of all
			sp_all.push_back(&(_species.back()));
		};

		// Create the interaction params
		Species *sp,*sp1,*sp2,*sp3;
		for (auto d: dims) {
			if (d.type==H) {

				// Create
				sp = _not_nullptr(_find_species(d.species1));
				_ixn_params.push_back(IxnParam(d.name,Hp,sp,d.guess));

				// Add to the species
				sp->set_h_ptr(&_ixn_params.back());

			} else if (d.type==J) { 

				// Create
				sp1 = _not_nullptr(_find_species(d.species1));
				sp2 = _not_nullptr(_find_species(d.species2));
				_ixn_params.push_back(IxnParam(d.name,Jp,sp1,sp2,d.guess));

				// Add to the species
				sp1->add_j_ptr(sp2,&_ixn_params.back());
				sp2->add_j_ptr(sp1,&_ixn_params.back());

			} else if (d.type==W) { 

				// Specific species or all?
				if (!d.all_species) {
					
					// Specific species
					
					// Create
					sp = _not_nullptr(_find_species(d.species1));
					_ixn_params.push_back(IxnParam(d.name,Wp,sp,d.guess));

					// Add to the species
					sp->add_w_ptr(&_ixn_params.back());

				} else {

					// All species

					// Create
					_ixn_params.push_back(IxnParam(d.name,Wp,sp_all,d.guess));

					// Add to all species
					for (auto spa: sp_all) {
						spa->add_w_ptr(&_ixn_params.back());
					};
				};

			} else if (d.type==B) {

				// Specific species or all?
				if (!d.all_species) {
					
					// Specific species

					// Create
					sp = _not_nullptr(_find_species(d.species1));
					_ixn_params.push_back(IxnParam(d.name,Bp,sp,d.guess));

					// No need to add to species, since it only affects activations

				} else {

					// All species

					// Create
					_ixn_params.push_back(IxnParam(d.name,Bp,sp_all,d.guess));

					// No need to add to species, since it only affects activations

				};
			} else if (d.type==K) {
				// Check that Lattice dim is 1 - only 1 is allowed!
				if (lattice_dim != 1) {
					std::cerr << "Error: triplets are currently only supported for lattice of dim 1" << std::endl;
					exit(EXIT_FAILURE);
				};

				// Create
				sp1 = _not_nullptr(_find_species(d.species1));
				sp2 = _not_nullptr(_find_species(d.species2));
				sp3 = _not_nullptr(_find_species(d.species3));
				_ixn_params.push_back(IxnParam(d.name,Kp,sp1,sp2,sp3,d.guess));

				// Add to species
				sp1->add_k_ptr(sp2,sp3,&_ixn_params.back());
				sp2->add_k_ptr(sp1,sp3,&_ixn_params.back());
				sp3->add_k_ptr(sp1,sp2,&_ixn_params.back());
			};
		};

		// Add a new counter for the dimensions
		IxnParam *ip;
		for (auto d: dims) {
			if (d.type==H) {
				add_counter(d.species1);
			} else if (d.type==J) {
				add_counter(d.species1,d.species2);
			} else if (d.type==K) {
				add_counter(d.species1,d.species2,d.species3);
			};
		};
	};
	BMLA::Impl::Impl(Impl&& other) : _latt(other._latt) {
		_copy(other);
		other._reset();
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
		_counters = other._counters;
		_latt = other._latt;
	};
	void BMLA::Impl::_reset() {
		_n_param = 0;
		_ixn_params.clear();
		_species.clear();
		_hidden_layer_exists = false;
		_hidden_units.clear();
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
		std::cout << "   Moments Initial: ";
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			std::cout << it->get_moment(IxnParam::AWAKE) << " ";
		};
		if (new_line) {
			std::cout << std::endl;
		};
		std::cout << "   Moments Final: ";
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			std::cout << it->get_moment(IxnParam::ASLEEP) << " ";
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

	void BMLA::Impl::add_counter(std::string s) {
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
	void BMLA::Impl::add_counter(std::string s1, std::string s2) {
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
	void BMLA::Impl::add_counter(std::string s1, std::string s2, std::string s3) {
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
	void BMLA::Impl::add_counter(std::string s1, std::string s2, std::string s3, std::string s4) {
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

	void BMLA::Impl::add_hidden_unit(std::vector<std::vector<int>> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params) {
		// Find sites indicated by connections
		std::vector<Site*> conns = _get_sites(lattice_idxs);
		// Make
		_add_hidden_unit(conns, w_params, b_params);
	};
	void BMLA::Impl::add_hidden_unit(std::vector<int> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params) {
		// Find sites indicated by connections
		std::vector<Site*> conns = _get_sites(lattice_idxs);
		// Make
		_add_hidden_unit(conns, w_params, b_params);
	};

	/********************
	Add hidden unit internal
	********************/

	void BMLA::Impl::_add_hidden_unit(std::vector<Site*> conn_sites, std::vector<std::string> w_params, std::vector<std::string> b_params) {

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

		// Combine conn_sites and the ixn params
		std::vector<std::pair<Site*,std::vector<IxnParam*>>> conns;
		for (auto c: conn_sites) {
			conns.push_back(std::make_pair(c,ip_w));
		};

		// Make hidden unit
		_hidden_units.push_back(HiddenUnit(conns,ip_b));

		// Go through lattice sites
		std::vector<Species*> sp_vec;
		for (auto c: conns) {
			// Go through ixn params W
			for (auto ipw: c.second) {
				// Get the species associated with this ixn param
				sp_vec = ipw->get_species();
				// Go through the species
				for (auto sp: sp_vec) {
					// Add to the site that this species has a conn to a hidden unit
					c.first->hidden_conns[sp].push_back(std::make_pair(&_hidden_units.back(),ip_w));
				};

				// Add to ixn param
				ipw->add_visible_hidden_connection(c.first,&_hidden_units.back());
			};
		};

		// Go through the biases b
		for (auto b: ip_b) {
			// Add this hidden unit
			b->add_hidden_unit(&_hidden_units.back());
		};
	};

	/********************
	Validate hidden layer
	********************/

	void BMLA::Impl::validate_hidden() const {
		for (auto it=_hidden_units.begin(); it != _hidden_units.end(); it++) {
			it->print_conns(true);
		};
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
			it->reset();
			it->moments_reset(IxnParam::AWAKE);
			it->moments_reset(IxnParam::ASLEEP);

			// Reset solution trajs as needed
			if (options.track_soln_traj) {
				it->set_track_soln_traj(true);
				it->reset_soln_traj();
			};
		};

		// Write the initial point for the solution
		if (options.write_soln_traj) {
			write(options.fname_write_soln_traj,true,0,false);
		};

		// Iterate over optimization steps
		for (int i_opt=0; i_opt<n_opt; i_opt++)
		{
			if (options.verbose) {
				std::cout << "Opt step: " << i_opt << " / " << n_opt << std::endl;
			};

			// Check MSE to see if quit
			if (options.mse_quit_mode) {
				if (_get_mse() < options.mse_quit) {
					// Quit!
					std::cout << "--- MSE is low enough: " << _get_mse() << " < " << options.mse_quit << " quitting! ---" << std::endl;
					break;
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
				if (options.awake_visible_are_binary) {
					// Binary
					_latt.read_from_file(fnames[randI(0,fnames.size()-1)],true);
				} else {
					// Probabilistic
					_latt.read_from_file(fnames[randI(0,fnames.size()-1)],false);
					// Sample -> binary
					// _latt.binarize();
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
				if (i_opt==0) {
					// Make new
					write_moments(options.fname_write_moment_traj,i_opt+1,false);
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

		std::list<Species>::iterator sp_it1,sp_it2,sp_it3;

		// Start by populating lattice randomly

		// Random number of initial particles (min is 1, max is box vol)
		int n = randI(1, pow(_latt.box_length(),_latt.dim()));

		// Random initial counts
		// Don't populate too much, else this is hard to find empty sites to place mols!
		// At most half the lattice is filled
		int n_possible = pow(_latt.box_length(),_latt.dim()) / 2;
		std::map<Species*,int> counts;
		for (std::list<Species>::iterator sp=_species.begin(); sp != _species.end(); sp++) {
			counts[&(*sp)] = randI(0,n_possible);
			n_possible -= counts[&(*sp)];
			if (n_possible < 0) { n_possible = 0; };
		};

		// Populate at random positions
		_latt.populate_randomly(counts);

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
			add_counter(sp);
			ctr = _not_nullptr(_find_ctr_by_species(sp));
			ctr->set_report_during_sampling(true);
		};
		for (auto sp_nbrs: options.report_nns) {
			_latt.init_nn_structure();
			add_counter(sp_nbrs.first,sp_nbrs.second);
			ctr = _not_nullptr(_find_ctr_by_species(sp_nbrs.first,sp_nbrs.second));
			ctr->set_report_during_sampling(true);
		};
		for (auto sp_triplets: options.report_triplets) {
			_latt.init_triplet_structure();
			add_counter(sp_triplets[0],sp_triplets[1],sp_triplets[2]);
			ctr = _not_nullptr(_find_ctr_by_species(sp_triplets[0],sp_triplets[1],sp_triplets[2]));
			ctr->set_report_during_sampling(true);
		};
		for (auto sp_quartics: options.report_quartics) {
			_latt.init_quartic_structure();
			add_counter(sp_quartics[0],sp_quartics[1],sp_quartics[2],sp_quartics[3]);
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
	IxnParam* BMLA::Impl::_find_ixn_param_by_species(IxnParamType type, std::string s) {
		if (type != Hp && type != Wp && type != Bp) {
			std::cerr << "ERROR: number of species does not match dimension type" << std::endl;
			exit(EXIT_FAILURE);
		};
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			if (it->is_type_with_species(type,s)) {
				return &*it;
			};
		};
		return nullptr;
	};
	IxnParam* BMLA::Impl::_find_ixn_param_by_species(IxnParamType type, std::string s1, std::string s2) {
		if (type != Jp) {
			std::cerr << "ERROR: number of species does not match dimension type" << std::endl;
			exit(EXIT_FAILURE);
		};
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			if (it->is_type_with_species(type,s1, s2)) {
				return &*it;
			};
		};
		return nullptr;
	};
	IxnParam* BMLA::Impl::_find_ixn_param_by_species(IxnParamType type, std::string s1, std::string s2, std::string s3) {
		if (type != Kp) {
			std::cerr << "ERROR: number of species does not match dimension type" << std::endl;
			exit(EXIT_FAILURE);
		};
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			if (it->is_type_with_species(type,s1,s2,s3)) {
				return &*it;
			};
		};
		return nullptr;
	};
	IxnParam* BMLA::Impl::_find_ixn_param_for_any_species(IxnParamType type) {
		if (type != Wp && type != Bp) {
			std::cerr << "ERROR: only Wp or Bp for any species." << std::endl;
			exit(EXIT_FAILURE);
		};
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			if (it->is_type_for_any_species(type)) {
				return &*it;
			};
		};
		return nullptr;
	};

	/****************************************
	BMLA IMPL forwards
	****************************************/

	// Constructor
	BMLA::BMLA(std::vector<Dim> dims, int box_length, int lattice_dim) : _impl(new Impl(dims,box_length,lattice_dim)) {};
	BMLA::BMLA(BMLA&& other) = default; // movable but no copies
    BMLA& BMLA::operator=(BMLA&& other) = default; // movable but no copies
	BMLA::~BMLA() = default;

	// Set a parameter for dim
	void BMLA::set_param_for_dim(std::string dim_name, double val) {
		_impl->set_param_for_dim(dim_name,val);
	};

	// Add hidden unit for any dim
	void BMLA::add_hidden_unit(std::vector<std::vector<int>> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params) {
		_impl->add_hidden_unit(lattice_idxs, w_params, b_params);
	};
	// 1D specific
	void BMLA::add_hidden_unit(std::vector<int> lattice_idxs, std::vector<std::string> w_params, std::vector<std::string> b_params) {
		_impl->add_hidden_unit(lattice_idxs,w_params,b_params);
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

	// Add a counter for some species or nns
	void BMLA::add_counter(std::string s) {
		_impl->add_counter(s);
	};
	void BMLA::add_counter(std::string s1, std::string s2) {
		_impl->add_counter(s1,s2);
	};
	void BMLA::add_counter(std::string s1, std::string s2, std::string s3) {
		_impl->add_counter(s1,s2,s3);
	};
	void BMLA::add_counter(std::string s1, std::string s2, std::string s3, std::string s4) {
		_impl->add_counter(s1,s2,s3,s4);
	};
};


