#include "dynamic_boltzmann.hpp"
#include "general.hpp"
#include <iostream>

using namespace DynamicBoltzmann;

int main() {

	// Dimensions vec
	std::vector<Dim> dims;
	dims.push_back(Dim("hA",DimType::H,"A",{"hA","jAA"},-4.0,0.0,20,-3.40145));
	dims.push_back(Dim("jAA",DimType::J,"A","A",{"hA","jAA"},0.0,4.0,20,2.74435));

	// Times
	double t_max = 1.0;
	int n_t = 101;

	// Opt params
	int batch_size = 10;
	int n_annealing = 10000;
	int box_length = 1000;
	double dopt = 0.01;
	int n_opt = 100;
	
	// Init
	std::cout << "Initializing..." << std::flush;
	OptProblem opt(dims, {"A"}, t_max, n_t, batch_size, n_annealing, box_length, dopt, n_opt, 1);
	std::cout << "ok." << std::endl;

	// Output
	opt.set_dir_io("bimol_annihilation_visible_data/");

	// Add filenames
	for (int i=1; i<=20; i++) {
		opt.add_fname("bimol_annihilation/lattice_v" + pad_str(i,2) + "/lattice/");
	};

	// Always use the same lattice for the batch
	bool use_same_lattice = true;

	// Solve
	std::cout << "Solving..." << std::endl;
	opt.solve(true, use_same_lattice);
	std::cout << "fin." << std::endl;

	return 0;
};