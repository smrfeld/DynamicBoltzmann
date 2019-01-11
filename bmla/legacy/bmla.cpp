#include "bmla.hpp"
#include <iostream>
#include <fstream>
#include "math.h"

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/********************
	Constructor
	********************/

	BMLA::BMLA(int box_length, int n_annealing, double dopt, std::vector<std::string> species, std::map<std::string,double> h0, std::map<std::string,std::map<std::string,double>> j0, int n_params, int n_opt, int n_batch) : _latt(box_length) {

		_n_annealing = n_annealing;
		_dopt = dopt;
		_n_params = n_params;
		_n_opt = n_opt;
		_n_batch = n_batch;

		// Make species
		for (auto s: species) {
			_species.push_back(Species(s));
			// Add to lattice
			_latt.add_species(&(_species.back()));
		};

		// Interaction params
		_t_fake=0; // always; fake time
		_soln_traj = new double*[1];
		_soln_traj[0] = new double[n_params];

		// Initialize moments
		_moms_init = new double[n_params];
		_moms_final = new double[n_params];

		// Initial interaction params
		for (auto m: h0) {
			_h0.push_back(m.second);
			// Find the species with this name
			for (auto its=_species.begin(); its!=_species.end(); its++) {
				if (its->name == m.first) {
					_h_species.push_back(&*its);
					// Link index
					its->_h_index = _h_species.size()-1;
					break;
				};
			};
		};
		for (auto m1: j0) {
			for (auto m2: m1.second) {
				_j0.push_back(m2.second);
				// Find the two species with these names
				bool done = false;
				for (auto its1=_species.begin(); its1!=_species.end(); its1++) {
					for (auto its2=_species.begin(); its2!=_species.end(); its2++) {
						if ( (its1->name == m1.first && its2->name == m2.first) || (its1->name == m2.first && its2->name == m1.first) ) {
							_j_species.push_back(std::make_pair(&*its1,&*its2));
							// Link index
							its1->_j_index[&*its2] = _h_species.size() + _j_species.size() - 1;
							its2->_j_index[&*its1] = _h_species.size() + _j_species.size() - 1;
							done = true;
							break;
						};
					};
					if (done) { break; };
				};
			};
		};

		// Link the species to their interaction params
		for (auto its=_species.begin(); its!=_species.end(); its++) {
			its->_soln_traj_ptr = &_soln_traj;
			its->_t_opt_ptr = &_t_fake;
		};
	};
	BMLA::~BMLA()
	{
		delete[] _soln_traj[0];
		delete[] _soln_traj;
		delete[] _moms_init;
		delete[] _moms_final;
	};

	/********************
	Solve for the h,j corresponding to a given lattice
	********************/

	void BMLA::solve(std::string fname, bool verbose)
	{
		// Initialize the params
		for (int i=0; i<_h0.size(); i++) {
			_soln_traj[0][i] = _h0[i];
		};
		for (int i=0; i<_j0.size(); i++) {
			_soln_traj[0][i+_h0.size()] = _j0[i];
		};

		// Declare
		std::pair<Species*,Species*> sp_pair;

		// Record the initial moments - these will never change!
		_latt.read_from_file(fname);
		for (int i=0; i<_n_params; i++) {
			if (i<_h_species.size())
			{
				_moms_init[i] = _h_species[i]->count;
			} else {
				sp_pair = _j_species[i-_h_species.size()];
				_moms_init[i] = sp_pair.first->nn_count[sp_pair.second];
			};
		};

		// Iterate over optimization steps
		for (int i_opt=0; i_opt<_n_opt; i_opt++)
		{
			if (verbose) {
				std::cout << "Opt step: " << i_opt << " / " << _n_opt << std::endl;
			};

			// Print out the current params
			if (verbose) {
				std::cout << "   Params: ";
				for (int i=0; i<_n_params; i++) {
					std::cout << _soln_traj[0][i] << " ";
				};
				std::cout << std::endl;
			};

			// Reset the final moments
			for (int i=0; i<_n_params; i++) {
				_moms_final[i] = 0.0;
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
				for (int i=0; i<_n_params; i++) {
					if (i<_h_species.size())
					{
						_moms_final[i] += _h_species[i]->count / _n_batch;
					} else {
						sp_pair = _j_species[i-_h_species.size()];
						_moms_final[i] += sp_pair.first->nn_count[sp_pair.second] / _n_batch;
					};
				};
			};
			if (verbose) {
				std::cout << std::endl;
			};

			// Print out the moments
			if (verbose) {
				std::cout << "   Moments Initial: ";
				for (int i=0; i<_n_params; i++) {
					std::cout << _moms_init[i] << " ";
				};
				std::cout << std::endl;
				std::cout << "   Moments Final: ";
				for (int i=0; i<_n_params; i++) {
					std::cout << _moms_final[i] << " ";
				};
				std::cout << std::endl;
				// MSE
				double mse=0.0;
				for (int i=0; i<_n_params; i++) {
					mse += abs(_moms_init[i]-_moms_final[i])/_moms_init[i];
				};
				std::cout << "   MSE: " << 100*mse/_n_params << "%" << std::endl;
			};

			// Update the params
			for (int i=0; i<_n_params; i++) {
				_soln_traj[0][i] += -1.0 * _dopt *(_moms_final[i] - _moms_init[i]);
			};

			// Compare final moments
			if (!verbose && i_opt == _n_opt-1) {
				std::cout << "   Moments Initial: ";
				for (int i=0; i<_n_params; i++) {
					std::cout << _moms_init[i] << " ";
				};
				std::cout << std::endl;
				std::cout << "   Moments Final: ";
				for (int i=0; i<_n_params; i++) {
					std::cout << _moms_final[i] << " ";
				};
				std::cout << std::endl;	
				// MSE
				double mse=0.0;
				for (int i=0; i<_n_params; i++) {
					mse += abs(_moms_init[i]-_moms_final[i])/_moms_init[i];
				};
				std::cout << "   MSE: " << 100*mse/_n_params << "%" << std::endl;
			};
		};
	};

	/********************
	Update the initial params
	********************/

	void BMLA::update_initial_params(std::map<std::string,double> h0, std::map<std::string,std::map<std::string,double>> j0) 
	{
		// Update m
		for (auto m: h0) {
			// Go through species
			for (int i=0; i<_h_species.size(); i++) {
				if (_h_species[i]->name == m.first) {
					// Update
					_h0[i] = m.second;
					break;
				};
			};
		};

		// Update j
		for (auto m1: j0) {
			for (auto m2: m1.second) {
				// Go through species
				for (int i=0; i<_j_species.size(); i++) {
					if ((_j_species[i].first->name == m1.first && _j_species[i].second->name == m2.first) || (_j_species[i].first->name == m2.first && _j_species[i].second->name == m1.first)) {
						// Update
						_j0[i] = m2.second;
						break;
					};
				};
			};
		};
	};

	/********************
	Write the solutions
	********************/

	void BMLA::write(std::string fname) {
		std::ofstream f;
		f.open (fname);
		std::pair<Species*,Species*> sp_pair;
		for (int i=0; i<_n_params; i++) {
			if (i<_h_species.size())
			{
				f << _h_species[i]->name << " " << _soln_traj[0][i] << "\n";
			} else {
				sp_pair = _j_species[i-_h_species.size()];
				f << sp_pair.first->name << " " << sp_pair.second->name << " " << _soln_traj[0][i] << "\n";
			};
		};
		f.close();	
	};

};


