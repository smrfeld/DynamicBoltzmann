#include "dynamic_boltzmann.hpp"
#include "general.hpp"
#include <iostream>

using namespace DynamicBoltzmann;

int main() {

	/********************
	Test optimization problem
	********************/

	// Number of dimensions
	int n_dim = 9;

	// Dimensions vec
	std::vector<Dim> dims;
	dims.push_back(Dim("hA",-2.0,0.2,30,H,"A"));
	dims.push_back(Dim("hB",-1.2,0.2,30,H,"B"));
	dims.push_back(Dim("hC",-1.5,0.6,30,H,"C"));
	dims.push_back(Dim("jAA",-0.2,2.0,30,J,"A","A"));
	dims.push_back(Dim("jAB",-1.0,1.0,30,J,"A","B"));
	dims.push_back(Dim("jAC",-1.0,0.3,30,J,"A","C"));
	dims.push_back(Dim("jBB",-0.5,2.0,30,J,"B","B"));
	dims.push_back(Dim("jBC",-1.0,1.0,30,J,"B","C"));
	dims.push_back(Dim("jCC",-0.2,1.5,30,J,"C","C"));

	// Fh
	dims[0].set_basis_func_dims(std::vector<std::string>{"hA","hB","hC"});
	dims[1].set_basis_func_dims(std::vector<std::string>{"hA","hB","hC"});
	dims[2].set_basis_func_dims(std::vector<std::string>{"hA","hB","hC"});

	// Jaa
	dims[3].set_basis_func_dims(std::vector<std::string>{"jAA","jAB","jAC"});
	// Jab
	dims[4].set_basis_func_dims(std::vector<std::string>{"jAA","jAB","jBB"});
	// Jac
	dims[5].set_basis_func_dims(std::vector<std::string>{"jAA","jCC","jAC"});
	// Jbb
	dims[6].set_basis_func_dims(std::vector<std::string>{"jAB","jBB","jBC"});
	// Jbc
	dims[7].set_basis_func_dims(std::vector<std::string>{"jBB","jBC","jCC"});
	// Jcc
	dims[8].set_basis_func_dims(std::vector<std::string>{"jBC","jAC","jCC"});

	// Initial conditions
	std::vector<double> init;
	init.push_back(-1.0); // hA
	init.push_back(-1.0); // hB
	init.push_back(-1.0); // hC
	init.push_back(0.0); // jAA
	init.push_back(0.0); // jAB
	init.push_back(0.0); // jAC
	init.push_back(0.0); // jBB
	init.push_back(0.0); // jBC
	init.push_back(0.0); // jCC

	// Times
	double t_max=1.0;
	int n_t = 100;

	// Opt params
	int batch_size = 5;
	int n_annealing = 3000;
	int box_length = 10;
	double dopt = 0.05;
	int n_opt = 100;
	
	// Init
	std::cout << "Initializing..." << std::flush;
	OptProblem opt(n_dim, t_max, n_t, batch_size, n_annealing, box_length, dopt, n_opt, dims, init, {"A","B","C"});
	std::cout << "ok." << std::endl;

	// Add filenames
	for (int i=0; i<100; i++) {
		opt.add_fname("rossler/lattice_v" + pad_str(i,2) + "/lattice/");
	};

	// Solve
	std::cout << "Solving..." << std::endl;
	opt.solve(true);
	std::cout << "fin." << std::endl;

	return 0;
};