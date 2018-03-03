#include "dynamic_boltzmann.hpp"
#include "general.hpp"
#include <iostream>

using namespace DynamicBoltzmann;

int main() {

	/********************
	Test optimization problem
	********************/

	// Number of dimensions
	int n_dim = 5;

	// Dimensions vec
	std::vector<Dim> dims;
	dims.push_back(Dim("hA",-3.0,0.15,16,H,"A"));
	dims.push_back(Dim("hB",-3.0,0.15,16,H,"B"));
	dims.push_back(Dim("jAA",-3.0,-0.05,16,J,"A","A"));
	dims.push_back(Dim("jAB",-3.0,-0.05,16,J,"A","B"));
	dims.push_back(Dim("jBB",-3.0,-0.05,16,J,"B","B"));

	dims[0].set_basis_func_dims(std::vector<std::string>{"hA","hB"});
	dims[1].set_basis_func_dims(std::vector<std::string>{"hA","hB"});
	// dims[0].set_basis_func_dims(std::vector<std::string>{"hA","hB","jAB"});
	// dims[1].set_basis_func_dims(std::vector<std::string>{"hA","hB","jAB"});
	dims[2].set_basis_func_dims(std::vector<std::string>{"jAA","jAB","jBB"});
	dims[3].set_basis_func_dims(std::vector<std::string>{"jAA","jAB","jBB"});
	dims[4].set_basis_func_dims(std::vector<std::string>{"jAA","jAB","jBB"});

	// Initial conditions
	std::vector<double> init;
	init.push_back(0.1); // hA
	init.push_back(0.1); // hB
	init.push_back(-0.1); // jAA
	init.push_back(-0.1); // jAB
	init.push_back(-0.1); // jBB

	// Times
	double t_max=1.0;
	int n_t = 101;

	// Opt params
	int batch_size = 10;
	int n_annealing = 300;
	int box_length = 10;
	double dopt = 0.5;
	int n_opt = 30;
	
	// Init
	std::cout << "Initializing..." << std::flush;
	OptProblem opt(n_dim, t_max, n_t, batch_size, n_annealing, box_length, dopt, n_opt, dims, init, {"A","B"});
	std::cout << "ok." << std::endl;

	// Add filenames
	for (int i=0; i<100; i++) {
		opt.add_fname("bimol_annihilation/lattice_v" + pad_str(i,2) + "/lattice/");
	};

	// Solve
	std::cout << "Solving..." << std::endl;
	opt.solve(true);
	std::cout << "fin." << std::endl;

	return 0;
};