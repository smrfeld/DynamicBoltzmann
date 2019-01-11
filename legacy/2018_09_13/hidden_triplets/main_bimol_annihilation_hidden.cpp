#include "../include/dynamic_boltzmann.hpp"
#include "../include/general.hpp"
#include <iostream>

using namespace DynamicBoltzmann;

int main() {

	// Dimensions vec
	std::vector<Dim> dims;
	dims.push_back(Dim("hA",DimType::H,"A",{"hA","wA"},-4.5,-2.5,20,-3.4886));
	dims.push_back(Dim("wA",DimType::W,"A",{"hA","wA"},-7.0,-4.0,20,-5.98862));

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
	opt.set_dir_io("bimol_annihilation_hidden_data/");

	// Add filenames
	for (int i=1; i<=20; i++) {
		opt.add_fname("bimol_annihilation/lattice_v" + pad_str(i,2) + "/lattice/");
	};

	// Setup the hidden layer
	// Connect every visible unit to a hidden unit
	// Each hidden unit is connected to 3 visible units, except at the ends where we only connect to 2
	std::vector<int> lattice_idxs;

	// End #1
	lattice_idxs.clear();
	lattice_idxs.push_back(1);
	lattice_idxs.push_back(2);
	opt.add_hidden_unit(lattice_idxs, "A");
	// End #2
	lattice_idxs.clear();
	lattice_idxs.push_back(box_length-1);
	lattice_idxs.push_back(box_length);
	opt.add_hidden_unit(lattice_idxs, "A");
	// Everything in the middle
	for (int idx = 2; idx<=box_length-1; idx++) {
		lattice_idxs.clear();
		lattice_idxs.push_back(idx-1);
		lattice_idxs.push_back(idx);
		lattice_idxs.push_back(idx+1);
		opt.add_hidden_unit(lattice_idxs, "A");
	};

	// Always use the same lattice for the batch
	bool use_same_lattice = true;

	// Solve
	std::cout << "Solving..." << std::endl;
	opt.solve(true, use_same_lattice);
	std::cout << "fin." << std::endl;

	return 0;
};