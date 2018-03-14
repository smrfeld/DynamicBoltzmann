#include "bmla.hpp"
#include <iostream>
#include <fstream>
#include "math.h"
#include <sstream>
#include <iomanip>

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Struct for specifying a dimension
	****************************************/

	Dim::Dim(std::string name, DimType type, std::string species, double guess) : Dim(name, type, species, "", guess) {};
	Dim::Dim(std::string name, DimType type, std::string species1, std::string species2, double guess)
	{
		if ((type==H && species2 != "") || (type==J && species2 == "")) {
			std::cerr << "ERROR! Dim specification is incorrect." << std::endl;
			exit(EXIT_FAILURE);
		};

		this->name = name;
		this->type = type;
		this->species1 = species1;
		this->species2 = species2;
		this->guess = guess;
	};

	/********************
	Constructor
	********************/

	BMLA::BMLA(std::vector<Dim> dims, std::vector<std::string> species, int batch_size, int n_annealing, int box_length, double dopt, int n_opt, int lattice_dim) : _latt(lattice_dim,box_length) {
		// Set parameters
		_n_param = dims.size();
		_dopt = dopt;
		_n_opt = n_opt;
		_n_annealing = n_annealing;
		_n_batch = batch_size;
		_mse_quit_mode = false;
		_mse_quit = 0.;
		_l2_reg = false;
		_lambda = 0.;

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
			} else { 
				_ixn_params.push_back(IxnParam(d.name,Jp,_find_species(d.species1),_find_species(d.species2),d.guess));
			};
		};

		// Add the interaction params and time ptr to the species
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
	};
	BMLA::BMLA(const BMLA& other) {
		_copy(other);
	};
	BMLA::BMLA(BMLA&& other) {
		_copy(other);
	};
	BMLA& BMLA::operator=(const BMLA& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
    BMLA& BMLA::operator=(BMLA&& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
	BMLA::~BMLA()
	{
		_clean_up();
	};

	void BMLA::_clean_up() {
		// Nothing...
	};
	void BMLA::_copy(const BMLA& other) {
		_n_param = other._n_param;
		_ixn_params = other._ixn_params;
		_species = other._species;
		_n_batch = other._n_batch;
		_n_annealing = other._n_annealing;
		_latt = other._latt;
		_dopt = other._dopt;
		_n_opt = other._n_opt;
		_mse_quit_mode = other._mse_quit_mode;
		_mse_quit = other._mse_quit;
		_l2_reg = other._l2_reg;
		_lambda = other._lambda;
	};
	void BMLA::_copy(BMLA&& other) {
		_n_param = other._n_param;
		_ixn_params = other._ixn_params;
		_species = other._species;
		_n_batch = other._n_batch;
		_n_annealing = other._n_annealing;
		_latt = other._latt;
		_dopt = other._dopt;
		_n_opt = other._n_opt;
		_mse_quit_mode = other._mse_quit_mode;
		_mse_quit = other._mse_quit;
		_l2_reg = other._l2_reg;
		_lambda = other._lambda;
		// Clear other
		other._n_param = 0;
		other._ixn_params.clear();
		other._species.clear();
		other._n_batch = 0;
		other._n_annealing = 0;
		other._latt = Lattice();
		other._dopt = 0.;
		other._n_opt = 0;
		other._mse_quit_mode = false;
		other._mse_quit = 0.;
		other._l2_reg = false;
		other._lambda = 0.;
	};

	/********************
	Print (compare) moments
	********************/

	void BMLA::_print_ixn_params(bool new_line) const {
		std::cout << "   Ixn params: ";
 		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			std::cout << it->get() << " ";
		};
		if (new_line) {
			std::cout << std::endl;
		};
	};

	void BMLA::_print_moments() const {
		std::cout << "   Moments Initial: ";
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			std::cout << it->get_moment(IxnParam::AWAKE) << " ";
		};
		std::cout << std::endl;
		std::cout << "   Moments Final: ";
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			std::cout << it->get_moment(IxnParam::ASLEEP) << " ";
		};
		std::cout << std::endl;
	};
	void BMLA::_print_mse(bool new_line) const {
		std::cout << "   MSE: " << _get_mse() << "%";
		if (new_line) {
			std::cout << std::endl;
		};
	};

	/********************
	Get MSE
	********************/

	double BMLA::_get_mse() const {
		double mse=0.0;
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			mse += abs(it->moments_diff())/it->get_moment(IxnParam::AWAKE);
		};
		return 100*mse/_n_param;
	};

	/********************
	Set and turn on l2 reg
	********************/

	void BMLA::set_l2_reg(double lambda) {
		_l2_reg = true;
		_lambda = lambda;
	};

	/********************
	Set and turn on MSE quit mode
	********************/

	void BMLA::set_mse_quit(double mse_quit) {
		_mse_quit_mode = true;
		_mse_quit = mse_quit;
	};

	/********************
	Solve for the h,j corresponding to a given lattice
	********************/

	void BMLA::solve(std::string fname, bool verbose)
	{
		// Reset the params to the guesses, and the moments to 0
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			it->reset();
			it->moments_reset(IxnParam::AWAKE);
			it->moments_reset(IxnParam::ASLEEP);
		};

		// Record the initial moments - these will never change!
		_latt.read_from_file(fname);
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			it->moments_retrieve(IxnParam::AWAKE);
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

			// Reset the asleep moments
			for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
				it->moments_reset(IxnParam::ASLEEP);
			};

			// Go over the batch
			if (verbose) { std::cout << "   " << std::flush; };
			for (int i_batch=0; i_batch<_n_batch; i_batch++)
			{
				if (verbose) {
					std::cout << "." << std::flush;
				};

				// Reset the lattice by reading it in
				_latt.read_from_file(fname);

				// Anneal
				_latt.anneal(_n_annealing);

				// Update the final moments
				for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
					it->moments_retrieve(IxnParam::ASLEEP, _n_batch);
				};
			};

			// Print out the MSE
			if (verbose) {
				_print_ixn_params(false);
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

	void BMLA::read(std::string fname) 
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

	void BMLA::write(std::string fname) {
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

	Species* BMLA::_find_species(std::string name) {
		for (auto it=_species.begin(); it!=_species.end(); it++) {
			if (it->name() == name) {
				return &*it;
			};
		};
		std::cerr << "ERROR: could not find species: " << name << std::endl;
		exit(EXIT_FAILURE);
	};
	IxnParam* BMLA::_find_ixn_param(std::string name) {
		for (auto it=_ixn_params.begin(); it!=_ixn_params.end(); it++) {
			if (it->name() == name) {
				return &*it;
			};
		};
		std::cerr << "ERROR: could not find ixn param: " << name << std::endl;
		exit(EXIT_FAILURE);
	};

};


