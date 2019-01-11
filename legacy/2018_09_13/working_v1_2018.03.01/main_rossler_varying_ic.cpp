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
	dims.push_back(Dim("hA",-2.0,2.0,21,H,"A"));
	dims.push_back(Dim("hB",-2.0,2.0,21,H,"B"));
	dims.push_back(Dim("hC",-2.0,2.0,21,H,"C"));
	dims.push_back(Dim("jAA",-2.0,2.0,21,J,"A","A"));
	dims.push_back(Dim("jAB",-2.0,2.0,21,J,"A","B"));
	dims.push_back(Dim("jAC",-2.0,2.0,21,J,"A","C"));
	dims.push_back(Dim("jBB",-2.0,2.0,21,J,"B","B"));
	dims.push_back(Dim("jBC",-2.0,2.0,21,J,"B","C"));
	dims.push_back(Dim("jCC",-2.0,2.0,21,J,"C","C"));

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

	// Times
	double t_max=0.7;
	int n_t = 70;

	// Opt params
	int batch_size = 10;
	int n_annealing = 1000000;
	int box_length = 10;
	double dopt = 1.0;
	int n_opt = 30;
	
	// Init
	std::cout << "Initializing..." << std::flush;
	OptProblem opt(n_dim, t_max, n_t, batch_size, n_annealing, box_length, dopt, n_opt, dims, std::vector<double>(), {"A","B","C"});
	std::cout << "ok." << std::endl;

	// Add filenames
	for (int i=0; i<142; i++) {
		opt.add_fname("rossler/" + pad_str(i,3) + "/");
	};

	// Solve
	std::cout << "Solving..." << std::endl;
	opt.solve_varying_ic(true);
	std::cout << "fin." << std::endl;

	return 0;
};