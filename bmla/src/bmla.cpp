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
		this->species1 = species1;
		this->species2 = species2;
		this->species3 = species3;
		this->guess = guess;
	};

	/****************************************
	BMLA - IMPLEMENTATION
	****************************************/

	// Declare
	struct MomRet;

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

		// Filename to write soln traj to
		std::string _fname_write_soln;

		// Filename to write moment traj to
		std::string _fname_write_moments;

		// Whether or not the visible layer is binary
		bool _visible_binary;

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
		IxnParam* _find_ixn_param_by_name(std::string name, bool enforce_success=true);
		IxnParam* _find_ixn_param_j_by_species(std::string species_name_1, std::string species_name_2, bool enforce_success=true);
		IxnParam* _find_ixn_param_w_by_species(std::string species_name, bool enforce_success=true);
		IxnParam* _find_ixn_param_k_by_species(std::string species_name_1, std::string species_name_2, std::string species_name_3, bool enforce_success=true);
		IxnParam* _find_ixn_param_b_by_species(std::string species_name, bool enforce_success=true);

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

		// Set a parameter for dim
		void set_param_for_dim(std::string dim_name, double val);

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

		// Set flag that we should write out the trajectory of parameters solved over the opt steps
		void set_write_soln_traj(std::string fname);

		// Set flag that we should write out the trajectory of moments solved for
		void set_write_moment_traj(std::string fname);

		// Set flag the visible units are binary
		// Default = true
		void set_binary_visible(bool flag);

		// Solve for the h,j corresponding to a given lattice
		void solve(std::vector<std::string> fnames, BinaryChoices binary_choices, bool verbose=false);

		// At the current ixns params, sample and report the specified moments
		// batch size = 1
		void sample(int n_cd_steps, bool write, std::string fname, BinaryChoices binary_choices, bool report_h, bool report_j, bool report_k, bool verbose);
		// given batch size
		void sample(int batch_size, int n_cd_steps, bool write, std::string fname, BinaryChoices binary_choices, bool report_h, bool report_j, bool report_k, bool verbose);
		// Internal
		MomRet _sample(int n_cd_steps, BinaryChoices binary_choices, bool report_h, bool report_j, bool report_k);

		// Update the initial params
		void read(std::string fname);

		// Write out the solutions
		void write(std::string fname, bool append, int opt_step);
		void write_moments(std::string fname, bool append, int idx1);
		void write_moments(std::string fname, bool append, int idx1, int idx2);
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
		_fname_write_soln = "";
		_fname_write_moments = "";
		_visible_binary = true;

		// Tell the lattice about what dims exist
		for (auto d: dims) {
			if (d.type==H) { 
				_latt.set_exists_h(true);
			} else if (d.type==J) {
				_latt.set_exists_j(true);
			} else if (d.type==K) {
				_latt.set_exists_k(true);
			} else if (d.type==W) {
				_hidden_layer_exists = true;
				_latt.set_exists_w(true);
			} else if (d.type==B) {
				_hidden_layer_exists = true;
				// No need to tell lattice, since it only affects hidden unit activation and not sampling
			};
		};

		// Create the species and add to the lattice
		for (auto s: species) {
			_species.push_back(Species(s));

			// Add to the lattice
			_latt.add_species(&(_species.back()));
		};

		// Tell the species about the other species to initialize counts
		for (auto itsp = _species.begin(); itsp!=_species.end(); itsp++) {
			itsp->init_counts(_species);
		};

		// Create the interaction params
		for (auto d: dims) {
			if (d.type==H) { 
				_ixn_params.push_back(IxnParam(d.name,Hp,_find_species(d.species1),d.guess));
			} else if (d.type==J) { 
				_ixn_params.push_back(IxnParam(d.name,Jp,_find_species(d.species1),_find_species(d.species2),d.guess));
			} else if (d.type==W) { 
				_ixn_params.push_back(IxnParam(d.name,Wp,_find_species(d.species1),d.guess));
			} else if (d.type==B) {
				_ixn_params.push_back(IxnParam(d.name,Bp,_find_species(d.species1),d.guess));
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
			} else if (d.type==K) {
				sp1 = _find_species(d.species1);
				sp2 = _find_species(d.species2);
				sp3 = _find_species(d.species3);
				sp1->add_k_ptr(sp2,sp3,ip_ptr);
				sp2->add_k_ptr(sp1,sp3,ip_ptr);
				sp3->add_k_ptr(sp1,sp2,ip_ptr);
			};
			// No need for Bp, since it only affects the hidden units (activation)
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
		_fname_write_soln = other._fname_write_soln;
		_fname_write_moments = other._fname_write_moments;
		_visible_binary = other._visible_binary;
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
		_fname_write_soln = "";
		_fname_write_moments = "";
		_visible_binary = true;
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
	Set a parameter for dim
	********************/

	void BMLA::Impl::set_param_for_dim(std::string dim_name, double val) {
		// Find
		IxnParam *ip = _find_ixn_param_by_name(dim_name);

		// Guaranteed not to be null
		ip->set_val(val);
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

		// See if a bias exists for hidden units with this species
		ip = _find_ixn_param_b_by_species(species,false); // can fail
		if (ip) {
			// Tell the bias that this hidden unit exists
			ip->add_hidden_unit(&_hidden_units.back());

			// Tell the hidden unit that this is it's bias
			_hidden_units.back().set_bias(ip);
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

	// Set flag that we should write out the trajectory of parameters solved over the opt steps
	void BMLA::Impl::set_write_soln_traj(std::string fname) {
		_fname_write_soln = fname;
	};

	// Set flag that we should write out the trajectory of moments solved for
	void BMLA::Impl::set_write_moment_traj(std::string fname) {
		_fname_write_moments = fname;
	};

	// Set flag the visible units are binary
	// Default = true
	void BMLA::Impl::set_binary_visible(bool flag) {
		_visible_binary = flag;
	};

	/********************
	Solve for the h,j corresponding to a given lattice
	********************/

	void BMLA::Impl::solve(std::vector<std::string> fnames, BinaryChoices binary_choices, bool verbose)
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

		// Write the initial point for the solution
		if (_fname_write_soln != "") {
			write(_fname_write_soln,false,0);
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
				if (_visible_binary) {
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
						// When using data, always use binary units
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
					if (cd_step != _n_cd_steps-1) {
						_latt.sample(binary_choices.asleep_visible_are_binary);
					} else {
						_latt.sample(binary_choices.asleep_final_visible_are_binary);
					};

					// Activate hidden
					if (_hidden_layer_exists) {
						for (auto ithu = _hidden_units.begin(); ithu != _hidden_units.end(); ithu++) {
							if (cd_step != _n_cd_steps-1) {
								ithu->activate(binary_choices.asleep_hidden_are_binary);
							} else {
								ithu->activate(binary_choices.asleep_final_hidden_are_binary);
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

			// Write the moments if needed
			if (_fname_write_moments != "") {
				if (i_opt==0) {
					// Make new
					write_moments(_fname_write_moments,false,i_opt+1);
				} else {
					// Append
					write_moments(_fname_write_moments,true,i_opt+1);
				};
			};

			// Update the params
			for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
				it->update(_dopt,_l2_reg,_lambda);
			};

			// Write the new solution (append)
			if (_fname_write_soln != "") {
				write(_fname_write_soln,true,i_opt+1);
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

	// Structure to return some moments
	struct MomRet {
		std::map<Species*,std::vector<double>> counts;
		std::map<Species*, std::map<Species*, std::vector<double>>> nns;
		std::map<Species*, std::map<Species*, std::map<Species*, std::vector<double>>>> triplets;

		// Incrementing
		// Caution: no check if keys exist!
		void increment(MomRet other, double batch_size=1.0) {
			for (auto pr: other.counts) {
				for (int i=0; i<pr.second.size(); i++) {
					counts[pr.first][i] += pr.second[i] / batch_size;
				};
			};
			for (auto pr1: other.nns) {
				for (auto pr2: pr1.second) {
					for (int i=0; i<pr2.second.size(); i++) {
						nns[pr1.first][pr2.first][i] += pr2.second[i] / batch_size;
					};
				};
			};
			for (auto pr1: other.triplets) {
				for (auto pr2: pr1.second) {
					for (auto pr3: pr2.second) {
						for (int i=0; i<pr3.second.size(); i++) {
							triplets[pr1.first][pr2.first][pr3.first][i] += pr3.second[i] / batch_size;
						};
					};
				};
			};
		};

		// Ave of int vector
		double _ave_of_vector(std::vector<double>& v) const {
			if (v.size() == 0) { return 0.; };
			double t=0.0;
			for (auto x: v) {
				t += x;
			};
			return t/v.size();
		};


		// Printing
		void report_h() {
			std::cout << "--- h moments ---" << std::endl;
			for (auto pr: counts) {
				std::cout << pr.first->name() << " ave: " << _ave_of_vector(pr.second) << " final: " << pr.second.back() << std::endl;
			};
		};
		void report_j() {
			std::cout << "--- j moments ---" << std::endl;
			for (auto pr1: nns) {
				for (auto pr2: pr1.second) {
					std::cout << pr1.first->name() << " " << pr2.first->name() << " ave: " << _ave_of_vector(pr2.second) << " final: " << pr2.second.back() << std::endl;
				};
			};
		};
		void report_k() {
			std::cout << "--- k moments ---" << std::endl;
			for (auto pr1: triplets) {
				for (auto pr2: pr1.second) {
					for (auto pr3: pr2.second) {
						std::cout << pr1.first->name() << " " << pr2.first->name() << " " << pr3.first->name() << " ave: " << _ave_of_vector(pr3.second) << " final: " << pr3.second.back() << std::endl;
					};
				};
			};
		};

		// Writing
		void write_h(std::ofstream &f) {
			for (auto pr: counts) {
				f << "h " <<  pr.first->name();
				for (auto x: pr.second) {
					f << " " << x;
				};
				f << std::endl;
			};
		};
		void write_j(std::ofstream &f) {
			for (auto pr1: nns) {
				for (auto pr2: pr1.second) {
					f << "J " <<  pr1.first->name() << " " << pr2.first->name();
					for (auto x: pr2.second) {
						f << " " << x;
					};
					f << std::endl;
				};
			};
		};
		void write_k(std::ofstream &f) {
			for (auto pr1: triplets) {
				for (auto pr2: pr1.second) {
					for (auto pr3: pr2.second) {
						f << "K " <<  pr1.first->name() << " " << pr2.first->name() << " " << pr3.first->name();
						for (auto x: pr3.second) {
							f << " " << x;
						};
						f << std::endl;
					};
				};
			};
		};
	};

	// Internal sampling function
	MomRet BMLA::Impl::_sample(int n_cd_steps, BinaryChoices binary_choices, bool report_h, bool report_j, bool report_k) {

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
				// When using real data, always use binary states
				ithu->activate(true);
			};
		};

		// Calculate and store initial moments
		MomRet moms;
		if (report_h) {
			for (sp_it1 = _species.begin(); sp_it1 != _species.end(); sp_it1++) {
				moms.counts[&(*sp_it1)].push_back(sp_it1->count());
			};
		};
		if (report_j) {
			for (sp_it1 = _species.begin(); sp_it1 != _species.end(); sp_it1++) {
				for (sp_it2 = sp_it1; sp_it2 != _species.end(); sp_it2++) {
					moms.nns[&(*sp_it1)][&(*sp_it2)].push_back(sp_it1->nn_count(&(*sp_it2)));
				};
			};
		};
		if (report_k) {
			for (sp_it1 = _species.begin(); sp_it1 != _species.end(); sp_it1++) {
				for (sp_it2 = sp_it1; sp_it2 != _species.end(); sp_it2++) {
					for (sp_it3 = sp_it2; sp_it3 != _species.end(); sp_it3++) {
						moms.triplets[&(*sp_it1)][&(*sp_it2)][&(*sp_it3)].push_back(sp_it1->triplet_count(&(*sp_it2),&(*sp_it3)));
					};
				};
			};
		};

		// CD
		for (int cd_step=0; cd_step<n_cd_steps; cd_step++) {

			// Sample visibles
			if (cd_step != n_cd_steps-1) {
				_latt.sample(binary_choices.asleep_visible_are_binary);
			} else {
				_latt.sample(binary_choices.asleep_final_visible_are_binary);
			};

			// Activate hidden
			if (_hidden_layer_exists) {
				for (auto ithu = _hidden_units.begin(); ithu != _hidden_units.end(); ithu++) {
					if (cd_step != n_cd_steps-1) {
						ithu->activate(binary_choices.asleep_hidden_are_binary);
					} else {
						ithu->activate(binary_choices.asleep_final_hidden_are_binary);
					}
				};
			};

			// Calculate the moments
			if (report_h) {
				for (sp_it1=_species.begin(); sp_it1 != _species.end(); sp_it1++) {
					moms.counts[&(*sp_it1)].push_back(sp_it1->count());
				};
			};
			if (report_j) {
				for (sp_it1=_species.begin(); sp_it1 != _species.end(); sp_it1++) {
					for (sp_it2=sp_it1; sp_it2 != _species.end(); sp_it2++) {
						moms.nns[&(*sp_it1)][&(*sp_it2)].push_back(sp_it1->nn_count(&(*sp_it2)));
					};
				};
			};
			if (report_k) {
				for (sp_it1=_species.begin(); sp_it1 != _species.end(); sp_it1++) {
					for (sp_it2=sp_it1; sp_it2 != _species.end(); sp_it2++) {
						for (sp_it3=sp_it2; sp_it3 != _species.end(); sp_it3++) {
							moms.triplets[&(*sp_it1)][&(*sp_it2)][&(*sp_it3)].push_back(sp_it1->triplet_count(&(*sp_it2),&(*sp_it3)));
						};
					};
				};
			};
		};

		// Return
		return moms;
	};

	// batch size = 1
	void BMLA::Impl::sample(int n_cd_steps, bool write, std::string fname, BinaryChoices binary_choices, bool report_h, bool report_j, bool report_k, bool verbose) {

		// Check writing
		if (write && fname == "") {
			std::cerr << "ERROR: provide a filename for writing!" << std::endl;
			exit(EXIT_FAILURE);
		};

		// Sample
		MomRet moms = _sample(n_cd_steps,binary_choices,report_h,report_j,report_k);

		// Report the moments
		if (verbose)
		{
			if (report_h) {
				moms.report_h();
			};
			if (report_j) {
				moms.report_j();
			};
			if (report_k) {
				moms.report_k();
			};
		};

		// Write if needed
		if (write) {
			std::ofstream f;
			f.open(fname);
			if (report_h) {
				moms.write_h(f);
			};
			if (report_j) {
				moms.write_j(f);
			};
			if (report_k) {
				moms.write_k(f);
			};
			f.close();
		};
	};

	// given batch size
	void BMLA::Impl::sample(int batch_size, int n_cd_steps, bool write, std::string fname, BinaryChoices binary_choices, bool report_h, bool report_j, bool report_k, bool verbose) {

		// Check writing
		if (write && fname == "") {
			std::cerr << "ERROR: provide a filename for writing!" << std::endl;
			exit(EXIT_FAILURE);
		};

		// Store average
		MomRet moms_ave, moms;

		// Go through batches
		for (int batch_no=0; batch_no<batch_size; batch_no++)
		{
			std::cout << "." << std::flush;

			moms = _sample(n_cd_steps,binary_choices,report_h,report_j,report_k);

			// Append
			if (batch_no == 0) {
				moms_ave = moms;
			} else {
				moms_ave.increment(moms,1.0*batch_size);
			};
		};
		std::cout << std::endl;

		// Report the moments
		if (verbose)
		{
			if (report_h) {
				moms_ave.report_h();
			};
			if (report_j) {
				moms_ave.report_j();
			};
			if (report_k) {
				moms_ave.report_k();
			};
		};

		// Write if needed
		if (write) {
			std::ofstream f;
			f.open(fname);
			if (report_h) {
				moms_ave.write_h(f);
			};
			if (report_j) {
				moms_ave.write_j(f);
			};
			if (report_k) {
				moms_ave.write_k(f);
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
		if (f.is_open()) { // make sure we found it
			while (getline(f,line)) {
				if (line == "") { continue; };
				iss = std::istringstream(line);
			    iss >> ixn_name;
			    iss >> guess;
		    	// Add
			    _find_ixn_param_by_name(ixn_name)->set_guess(atof(guess.c_str()));
		    	ixn_name=""; guess="";
			};
		};
		f.close();
	};

	/********************
	Write the solutions
	********************/

	void BMLA::Impl::write(std::string fname, bool append, int opt_step) {
		std::ofstream f;
		if (append) {
			f.open(fname, std::ofstream::out | std::ofstream::app);
		} else {
			f.open(fname);
		};
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			if (opt_step != -1) {
				f << opt_step << " " << it->name() << " " << it->get() << "\n";
			} else {
				f << it->name() << " " << it->get() << "\n";
			};
		};
		f.close();	
	};

	void BMLA::Impl::write_moments(std::string fname, bool append, int idx1) {
		write_moments(fname,append,idx1,-1);
	};
	void BMLA::Impl::write_moments(std::string fname, bool append, int idx1, int idx2) {
		std::ofstream f;
		if (append) {
			f.open(fname, std::ofstream::out | std::ofstream::app);
		} else {
			f.open(fname);
		};
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			if (idx1 != -1 && idx2 != -1) {
				f << idx1 << " " << idx2 << " " << it->name() << " " << it->get_moment(IxnParam::AWAKE) << " " << it->get_moment(IxnParam::ASLEEP) << "\n";
			} else if (idx2 != -1) {
				f << idx2 << " " << it->name() << " " << it->get_moment(IxnParam::AWAKE) << " " << it->get_moment(IxnParam::ASLEEP) << "\n";
			} else if (idx1 != -1) {
				f << idx1 << " " << it->name() << " " << it->get_moment(IxnParam::AWAKE) << " " << it->get_moment(IxnParam::ASLEEP) << "\n";
			} else {
				f << it->name() << " " << it->get_moment(IxnParam::AWAKE) << " " << it->get_moment(IxnParam::ASLEEP) << "\n";
			};
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
	IxnParam* BMLA::Impl::_find_ixn_param_by_name(std::string name, bool enforce_success) {
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
	IxnParam* BMLA::Impl::_find_ixn_param_b_by_species(std::string species_name, bool enforce_success) {
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			if (it->is_b_with_species(species_name)) {
				return &*it;
			};
		};
		if (enforce_success) {
			std::cerr << "ERROR: could not find B ixn param for species: " << species_name << std::endl;
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

	// Set a parameter for dim
	void BMLA::set_param_for_dim(std::string dim_name, double val) {
		_impl->set_param_for_dim(dim_name,val);
	};

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

	// Set flag that we should write out the trajectory of parameters solved over the opt steps
	void BMLA::set_write_soln_traj(std::string fname) {
		_impl->set_write_soln_traj(fname);
	};

	// Set flag that we should write out the trajectory of moments solved for
	void BMLA::set_write_moment_traj(std::string fname) {
		_impl->set_write_moment_traj(fname);
	};

	// Set flag the visible units are binary
	// Default = true
	void BMLA::set_binary_visible(bool flag) {
		_impl->set_binary_visible(flag);
	};

	// Solve for the h,j corresponding to a given lattice
	void BMLA::solve(std::string fname, BinaryChoices binary_choices, bool verbose) {
		std::vector<std::string> fnames;
		fnames.push_back(fname);
		_impl->solve(fnames,binary_choices,verbose);
	};

	// Solve for the h,j corresponding to a given lattice
	void BMLA::solve(std::vector<std::string> fnames, BinaryChoices binary_choices, bool verbose) {
		_impl->solve(fnames,binary_choices,verbose);
	};

	// At the current ixns params, sample and report the specified moments
	// batch size = 1
	void BMLA::sample(int n_cd_steps, bool write, std::string fname, BinaryChoices binary_choices, bool report_h, bool report_j, bool report_k, bool verbose) {
		_impl->sample(n_cd_steps, write, fname, binary_choices, report_h,report_j,report_k,verbose);
	};
	// given batch size
	void BMLA::sample(int batch_size, int n_cd_steps, bool write, std::string fname, BinaryChoices binary_choices, bool report_h, bool report_j, bool report_k, bool verbose) {
		_impl->sample(batch_size, n_cd_steps, write, fname, binary_choices, report_h,report_j,report_k,verbose);
	};

	// Update the initial params
	void BMLA::read(std::string fname) {
		_impl->read(fname);
	};

	// Write out the solutions
	void BMLA::write(std::string fname, bool append, int opt_step) {
		_impl->write(fname,append,opt_step);
	};
};


