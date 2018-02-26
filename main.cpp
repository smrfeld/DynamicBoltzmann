#include "dynamic_boltzmann.hpp"
#include "general.hpp"
#include <iostream>

using namespace DynamicBoltzmann;

int main() {

	/********************
	Test optimization problem
	********************/

	// Number of dimensions
	int n_dim = 2;

	// Dimensions vec
	std::vector<Dim> dims;
	dims.push_back(Dim("h",-0.6,0.55,21,H,"A"));
	dims.push_back(Dim("j",-0.7,0.15,21,J,"A","A"));

	// Initial conditions
	std::vector<double> init;
	init.push_back(0.5); // h
	init.push_back(0.1); // j

	// Times
	double t_max=1.0;
	int n_t = 101;

	// Opt params
	int batch_size = 10;
	int n_annealing = 500;
	int box_length = 10;
	double dopt = 0.005;
	int n_opt = 100;
	
	// Init
	std::cout << "Initializing..." << std::flush;
	OptProblem opt(n_dim, t_max, n_t, batch_size, n_annealing, box_length, dopt, n_opt, dims, init, {"A"});
	std::cout << "ok." << std::endl;

	// Add filenames
	for (int i=0; i<100; i++) {
		opt.add_fname("annihilation/lattice_v" + pad_str(i,2) + "/lattice/");
	};

	// Solve
	std::cout << "Solving..." << std::endl;
	opt.solve(true);
	std::cout << "fin." << std::endl;

	return 0;
};