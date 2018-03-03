#include "dynamic_boltzmann.hpp"
#include "general.hpp"
#include <iostream>

using namespace DynamicBoltzmann;

int main() {

	/********************
	Test optimization problem
	********************/

	// Dimensions vec
	std::vector<Dim> dims;
	dims.push_back(Dim("hA",DimType::H,"A",{"hA","hB"},-3.0,0.15,16,0.1));
	dims.push_back(Dim("hB",DimType::H,"B",{"hA","hB"},-3.0,0.15,16,0.1));
	dims.push_back(Dim("jAA",DimType::J,"A","A",{"jAA","jAB","jBB"},-3.0,-0.05,16,-0.1));
	dims.push_back(Dim("jAB",DimType::J,"A","B",{"jAA","jAB","jBB"},-3.0,-0.05,16,-0.1));
	dims.push_back(Dim("jBB",DimType::J,"B","B",{"jAA","jAB","jBB"},-3.0,-0.05,16,-0.1));

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
	OptProblem opt(dims, {"A","B"}, t_max, n_t, batch_size, n_annealing, box_length, dopt, n_opt);
	std::cout << "ok." << std::endl;

	// Add filenames
	for (int i=0; i<100; i++) {
		opt.add_fname("bimol_annihilation/lattice_v" + pad_str(i,2) + "/lattice/");
	};

	// Validate
	opt.validate_setup();

	// Solve
	std::cout << "Solving..." << std::endl;
	opt.solve(true);
	std::cout << "fin." << std::endl;

	return 0;
};