#include <iostream>
#include "../../bmla/include/bmla.hpp"
#include "../../bmla/include/general.hpp"
#include <fstream>

using namespace DynamicBoltzmann;

int main() {

	int box_length = 1000;
	int n_cd_steps = 10;
	double dopt = 0.002;
	int n_opt = 1000;
	int n_batch = 1000;

	std::vector<Dim> dims;
	dims.push_back(Dim("h",DimType::H,"A",0.83461));
	dims.push_back(Dim("W",DimType::W,"A",-7.175));

	BMLA bmla(dims,{"A"},n_batch,box_length,dopt,n_opt,1);

	// Set the number of CD steps
	bmla.set_n_cd_steps(n_cd_steps);

	// Quit if the MSE goes below 1%
	bmla.set_mse_quit(0.1);

	// Turn on L2 reg
	// bmla.set_l2_reg(0.01);

	// Add hidden units

	// Connect every visible unit to a hidden unit
	// Each hidden unit is connected to 3 visible units, except at the ends where we only connect to 2
	std::vector<int> lattice_idxs;

	// End #1
	lattice_idxs.clear();
	lattice_idxs.push_back(1);
	lattice_idxs.push_back(2);
	bmla.add_hidden_unit(lattice_idxs, "A");
	// End #2
	lattice_idxs.clear();
	lattice_idxs.push_back(box_length-1);
	lattice_idxs.push_back(box_length);
	bmla.add_hidden_unit(lattice_idxs, "A");
	// Everything in the middle
	for (int idx = 2; idx<=box_length-1; idx++) {
		lattice_idxs.clear();
		lattice_idxs.push_back(idx-1);
		lattice_idxs.push_back(idx);
		lattice_idxs.push_back(idx+1);
		bmla.add_hidden_unit(lattice_idxs, "A");
	};

	// Iterate over the things
	int n_files = 10;
	for (int i=10; i<=n_files; i++) {
		std::cout << "--- File: " << i << " / " << n_files << " ---" << std::endl; 

		// Read the existing as a guess
		// bmla.read("../bimol_annihilation/lattice_v"+pad_str(i,2)+"/lattice/init.txt");		

		// Solve
		bmla.solve("../bimol_annihilation/lattice_v"+pad_str(i,2)+"/lattice/0000.txt",true);
		
		// Write
		bmla.write("../bimol_annihilation/lattice_v"+pad_str(i,2)+"/lattice/init_hidden.txt");
	};

	return 0;
};