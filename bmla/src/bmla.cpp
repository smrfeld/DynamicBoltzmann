#include "../include/bmla.hpp"
#include <iostream>
#include <fstream>
#include "math.h"
#include <sstream>
#include <iomanip>

#include "ixn_param.hpp"
#include "../include/general.hpp"

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Struct for specifying a dimension
	****************************************/

	Dim::Dim(std::string name, DimType type, std::string species, double guess) : Dim(name, type, species, "", "", guess) {};
	Dim::Dim(std::string name, DimType type, std::string species1, std::string species2, double guess) : Dim(name, type, species1, species2, "", guess) {};
	Dim::Dim(std::string name, DimType type, std::string species1, std::string species2, std::string species3, double guess)
	{
		if (type == H || type == W) {
			if (species1 == "" || species2 != "" || species3 != "") {
				std::cerr << "ERROR! Dim specification is incorrect for H or W." << std::endl;
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

		// Batch size
		int _n_batch;

		// Number of CD steps
		int _n_cd_steps;

		// The lattice to learn
		Lattice _latt;

		// Optimization step size
		double _dopt;

		// No opt steps
		int _n_opt;

		// If the MSE dips below this, quit
		bool _mse_quit_mode;
		double _mse_quit; // in percent

		// L2 reg
		bool _l2_reg;
		double _lambda;

		// Use a single lattice, irregardless of batch size
		bool _use_single_lattice;

		// Print
		void _print_ixn_params(bool new_line=true) const;
		void _print_moments(bool new_line=true) const;
		void _print_mse(bool new_line=true) const;

		// Get the mse
		double _get_mse() const;

		// Add a hidden unit
		void _add_hidden_unit(std::vector<Site*> conns, std::string species);

		// Search functions
		Species* _find_species(std::string name, bool enforce_success=true);
		IxnParam* _find_ixn_param(std::string name, bool enforce_success=true);
		IxnParam* _find_ixn_param_j_by_species(std::string species_name_1, std::string species_name_2, bool enforce_success=true);
		IxnParam* _find_ixn_param_w_by_species(std::string species_name, bool enforce_success=true);
		IxnParam* _find_ixn_param_k_by_species(std::string species_name_1, std::string species_name_2, std::string species_name_3, bool enforce_success=true);

		// Constructor helpers
		void _clean_up();
		void _copy(const Impl& other);
		void _reset();

	public:

		// Constructor
		Impl(std::vector<Dim> dims, std::vector<std::string> species, int batch_size, int box_length, double dopt, int n_opt, int lattice_dim=3);
		Impl(Impl&& other);
	    Impl& operator=(Impl&& other);
		~Impl();

		// Any dim
		void add_hidden_unit(std::vector<std::vector<int>> lattice_idxs, std::string species);
		// 1D specific
		void add_hidden_unit(std::vector<int> lattice_idxs, std::string species);

		// Set the number of CD steps (default = 1)
		void set_n_cd_steps(int n_steps);

		// Set and turn on l2 regularizer
		void set_l2_reg(double lambda);

		// Set and turn on MSE quit mode
		void set_mse_quit(double mse_quit);

		// Use a single lattice for training, irregardless of batch size
		void set_use_single_lattice(bool flag);

		// Solve for the h,j corresponding to a given lattice
		void solve(std::vector<std::string> fnames, bool verbose=false);

		// Update the initial params
		void read(std::string fname);

		// Write out the solutions
		void write(std::string fname);
	};

	/****************************************
	BMLA - IMPLEMENTATION DEFINITIONS
	****************************************/

	/********************
	Constructor
	********************/

	BMLA::Impl::Impl(std::vector<Dim> dims, std::vector<std::string> species, int batch_size, int box_length, double dopt, int n_opt, int lattice_dim) : _latt(lattice_dim,box_length) {
		// Set parameters
		_n_param = dims.size();
		_dopt = dopt;
		_n_opt = n_opt;
		_n_batch = batch_size;
		_mse_quit_mode = false;
		_mse_quit = 0.;
		_n_cd_steps = 1; // default
		_l2_reg = false;
		_lambda = 0.;
		_hidden_layer_exists = false;
		_use_single_lattice = false; // default

		// Create the species and add to the lattice
		for (auto s: species) {
			_species.push_back(Species(s));

			// Add to the lattice
			_latt.add_species(&(_species.back()));
		};

		// Create the interaction params
		for (auto d: dims) {
			if (d.type==H) { 
				_ixn_params.push_back(IxnParam(d.name,Hp,_find_species(d.species1),d.guess));
			} else if (d.type==J) { 
				_ixn_params.push_back(IxnParam(d.name,Jp,_find_species(d.species1),_find_species(d.species2),d.guess));
			} else if (d.type==W) { 
				_ixn_params.push_back(IxnParam(d.name,Wp,_find_species(d.species1),d.guess));
			} else if (d.type==K) {
				// Check that Lattice dim is 1 - only 1 is allowed!
				if (lattice_dim != 1) {
					std::cerr << "Error: triplets are currently only supported for lattice of dim 1" << std::endl;
					exit(EXIT_FAILURE);
				};

				_ixn_params.push_back(IxnParam(d.name,Kp,_find_species(d.species1),_find_species(d.species2),_find_species(d.species3),d.guess));
			};
		};

		// Add the interaction params and time ptr to the species
		Species *sp1=nullptr, *sp2=nullptr, *sp3=nullptr;
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
			} else if (d.type==W) {
				sp1 = _find_species(d.species1);
				sp1->set_w_ptr(ip_ptr);		
			} else if (d.type==K) {
				sp1 = _find_species(d.species1);
				sp2 = _find_species(d.species2);
				sp3 = _find_species(d.species3);
				sp1->add_k_ptr(sp2,sp3,ip_ptr);
				sp2->add_k_ptr(sp1,sp3,ip_ptr);
				sp3->add_k_ptr(sp1,sp2,ip_ptr);
			};
		};
		// Ensure the J,K of the species are complete - if there is no ixn param that descripes the coupling, add a nullptr entry in the dictionary - later check if nullptr, then return 0
		// I think this is faster - otherwise there would be no reason to do it
		for (auto itsp1 = _species.begin(); itsp1!=_species.end(); itsp1++) {
			for (auto itsp2 = _species.begin(); itsp2!=_species.end(); itsp2++) {
				// J
				ip_ptr = _find_ixn_param_j_by_species(itsp1->name(), itsp2->name(), false);
				if (!ip_ptr) {
					// It's null; add to both
					itsp1->add_j_ptr(&(*itsp2),nullptr);
					itsp2->add_j_ptr(&(*itsp1),nullptr);
				};

				// K
				for (auto itsp3 = _species.begin(); itsp3!=_species.end(); itsp3++) {
					ip_ptr = _find_ixn_param_k_by_species(itsp1->name(), itsp2->name(), itsp3->name(), false);
					if (!ip_ptr) {
						// It's null; add to all 3
						itsp1->add_k_ptr(&(*itsp2),&(*itsp3),nullptr);
						itsp2->add_k_ptr(&(*itsp1),&(*itsp3),nullptr);
						itsp3->add_k_ptr(&(*itsp1),&(*itsp2),nullptr);
					};
				};
			};
		};
	};
	BMLA::Impl::Impl(Impl&& other) {
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
		_n_batch = other._n_batch;
		_latt = other._latt;
		_dopt = other._dopt;
		_n_opt = other._n_opt;
		_mse_quit_mode = other._mse_quit_mode;
		_mse_quit = other._mse_quit;
		_n_cd_steps = other._n_cd_steps;
		_l2_reg = other._l2_reg;
		_lambda = other._lambda;
		_use_single_lattice = other._use_single_lattice;
	};
	void BMLA::Impl::_reset() {
		_n_param = 0;
		_ixn_params.clear();
		_species.clear();
		_hidden_layer_exists = false;
		_hidden_units.clear();
		_n_batch = 0;
		_latt = Lattice();
		_dopt = 0.;
		_n_opt = 0;
		_mse_quit_mode = false;
		_mse_quit = 0.;
		_n_cd_steps = 1;
		_l2_reg = false;
		_lambda = 0.;
		_use_single_lattice = false;
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
			mse += abs(it->moments_diff())/it->get_moment(IxnParam::AWAKE);
		};
		return 100*mse/_n_param;
	};

	/********************
	Set properties
	********************/

	void BMLA::Impl::add_hidden_unit(std::vector<std::vector<int>> lattice_idxs, std::string species) 
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
	void BMLA::Impl::add_hidden_unit(std::vector<int> lattice_idxs, std::string species) 
	{
		// Find sites indicated by connections
		std::vector<Site*> conns;
		for (auto c: lattice_idxs) {
			conns.push_back(_latt.get_site(c));
		};

		// Add
		_add_hidden_unit(conns,species);
	};
	void BMLA::Impl::_add_hidden_unit(std::vector<Site*> conns, std::string species)
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
		IxnParam *ip = _find_ixn_param_w_by_species(species);
		for (auto sptr: conns) {
			ip->add_visible_hidden_connection(sptr,&_hidden_units.back());
		};
	};

	/********************
	Set the number of CD steps (default = 1)
	********************/

	void BMLA::Impl::set_n_cd_steps(int n_steps) {
		_n_cd_steps = n_steps;
	};

	/********************
	Set and turn on l2 reg
	********************/

	void BMLA::Impl::set_l2_reg(double lambda) {
		_l2_reg = true;
		_lambda = lambda;
	};

	/********************
	Set and turn on MSE quit mode
	********************/

	void BMLA::Impl::set_mse_quit(double mse_quit) {
		_mse_quit_mode = true;
		_mse_quit = mse_quit;
	};

	/********************
	Use a single lattice for training, irregardless of batch size
	********************/

	void BMLA::Impl::set_use_single_lattice(bool flag) {
		_use_single_lattice = flag;
	};

	/********************
	Solve for the h,j corresponding to a given lattice
	********************/

	void BMLA::Impl::solve(std::vector<std::string> fnames, bool verbose)
	{
		// Check size of filenames
		if (_use_single_lattice == true && fnames.size() != 1) {
			std::cerr << "Error! In _use_single_lattice mode, only provide one filename!" << std::endl;
			exit(EXIT_FAILURE); 
		};

		// Reset the params to the guesses, and the moments to 0
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			it->reset();
			it->moments_reset(IxnParam::AWAKE);
			it->moments_reset(IxnParam::ASLEEP);
		};

		// Iterate over optimization steps
		for (int i_opt=0; i_opt<_n_opt; i_opt++)
		{
			if (verbose) {
				std::cout << "Opt step: " << i_opt << " / " << _n_opt << std::endl;
			};

			// Check MSE to see if quit
			if (_mse_quit_mode) {
				if (_get_mse() < _mse_quit) {
					// Quit!
					std::cout << "--- MSE is low enough: " << _get_mse() << " < " << _mse_quit << " quitting! ---" << std::endl;
					break;
				};
			};

			// Reset all moments
			for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
				it->moments_reset(IxnParam::AWAKE);
				it->moments_reset(IxnParam::ASLEEP);
			};

			// Go over the batch
			if (verbose) { std::cout << "   " << std::flush; };
			for (int i_batch=0; i_batch<_n_batch; i_batch++)
			{
				if (verbose) {
					std::cout << "." << std::flush;
				};

				// Reset the lattice by reading in a random
				_latt.read_from_file(fnames[randI(0,fnames.size()-1)]);
		 		
				// Activate hidden
				if (_hidden_layer_exists) {
					for (auto ithu = _hidden_units.begin(); ithu != _hidden_units.end(); ithu++) {
						// When using real data, always use binary states
						ithu->activate(true);
					};
				};

				// Record awake moments
				for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
					it->moments_retrieve(IxnParam::AWAKE, _n_batch);
				};

				// Do the contrastive divergence
				for (int cd_step=0; cd_step<_n_cd_steps; cd_step++)
				{
					// Sample
					_latt.sample(false); // probabilistic

					// Activate hidden
					if (_hidden_layer_exists) {
						for (auto ithu = _hidden_units.begin(); ithu != _hidden_units.end(); ithu++) {
							// When using reconstructions, always use raw probabilities
							ithu->activate(false);
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
					it->moments_retrieve(IxnParam::ASLEEP, _n_batch);
				};
				//std::cout << "------------------------------" << std::endl;

			};

			// Print out the MSE
			if (verbose) {
				_print_ixn_params(false);
				_print_moments(false);
				_print_mse();
			};

			// Update the params
			for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
				it->update(_dopt,_l2_reg,_lambda);
			};
		};

		// Report final
		std::cout << "--- Final ---" << std::endl;
		_print_ixn_params();
		_print_moments();
		_print_mse();
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
		if (f.is_open()) { // make sure we found it
			while (getline(f,line)) {
				if (line == "") { continue; };
				iss = std::istringstream(line);
			    iss >> ixn_name;
			    iss >> guess;
		    	// Add
			    _find_ixn_param(ixn_name)->set_guess(atof(guess.c_str()));
		    	ixn_name=""; guess="";
			};
		};
		f.close();
	};

	/********************
	Write the solutions
	********************/

	void BMLA::Impl::write(std::string fname) {
		std::ofstream f;
		f.open (fname);
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			f << it->name() << " " << it->get() << "\n";
		};
		f.close();	
	};

	/********************
	Search functions
	********************/

	Species* BMLA::Impl::_find_species(std::string name, bool enforce_success) {
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
	IxnParam* BMLA::Impl::_find_ixn_param(std::string name, bool enforce_success) {
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
	IxnParam* BMLA::Impl::_find_ixn_param_j_by_species(std::string species_name_1, std::string species_name_2, bool enforce_success) {
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
	IxnParam* BMLA::Impl::_find_ixn_param_w_by_species(std::string species_name, bool enforce_success) {
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
	IxnParam* BMLA::Impl::_find_ixn_param_k_by_species(std::string species_name_1, std::string species_name_2, std::string species_name_3, bool enforce_success) {
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			if (it->is_k_with_species(species_name_1,species_name_2,species_name_3)) {
				return &*it;
			};
		};
		if (enforce_success) {
			std::cerr << "ERROR: could not find K ixn param for species: " << species_name_1 << " " << species_name_2 << " " << species_name_3 << std::endl;
			exit(EXIT_FAILURE);
		} else {
			return nullptr;
		};
	};

	/****************************************
	BMLA IMPL forwards
	****************************************/

	// Constructor
	BMLA::BMLA(std::vector<Dim> dims, std::vector<std::string> species, int batch_size, int box_length, double dopt, int n_opt, int lattice_dim) : _impl(new Impl(dims,species,batch_size,box_length,dopt,n_opt,lattice_dim)) {};
	BMLA::BMLA(BMLA&& other) = default; // movable but no copies
    BMLA& BMLA::operator=(BMLA&& other) = default; // movable but no copies
	BMLA::~BMLA() = default;

	// Any dim
	void BMLA::add_hidden_unit(std::vector<std::vector<int>> lattice_idxs, std::string species) {
		_impl->add_hidden_unit(lattice_idxs,species);
	};
	// 1D specific
	void BMLA::add_hidden_unit(std::vector<int> lattice_idxs, std::string species) {
		_impl->add_hidden_unit(lattice_idxs,species);
	};

	// Set the number of CD steps (default = 1)
	void BMLA::set_n_cd_steps(int n_steps) {
		_impl->set_n_cd_steps(n_steps);
	};

	// Set and turn on l2 regularizer
	void BMLA::set_l2_reg(double lambda) {
		_impl->set_l2_reg(lambda);
	};

	// Set and turn on MSE quit mode
	void BMLA::set_mse_quit(double mse_quit) {
		_impl->set_mse_quit(mse_quit);
	};

	// Use a single lattice for training, irregardless of batch size
	void BMLA::set_use_single_lattice(bool flag) {
		_impl->set_use_single_lattice(flag);
	};

	// Solve for the h,j corresponding to a given lattice
	void BMLA::solve(std::string fname, bool verbose) {
		std::vector<std::string> fnames;
		fnames.push_back(fname);
		_impl->solve(fnames,verbose);
	};

	// Solve for the h,j corresponding to a given lattice
	void BMLA::solve(std::vector<std::string> fnames, bool verbose) {
		_impl->solve(fnames,verbose);
	};

	// Update the initial params
	void BMLA::read(std::string fname) {
		_impl->read(fname);
	};

	// Write out the solutions
	void BMLA::write(std::string fname) {
		_impl->write(fname);
	};
};


