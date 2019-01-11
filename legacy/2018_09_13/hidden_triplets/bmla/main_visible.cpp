#include <iostream>
#include "../../bmla/include/bmla.hpp"
#include "../../bmla/include/general.hpp"
#include <fstream>

using namespace DynamicBoltzmann;

int main() {

	int box_length = 1000;
	int n_annealing = 100000;
	double dopt = 0.002;
	int n_opt = 100;
	int n_batch = 20;

	std::vector<Dim> dims;
	dims.push_back(Dim("h",DimType::H,"A",-3.4072));
	dims.push_back(Dim("J",DimType::J,"A","A",2.7214));

	BMLA bmla(dims,{"A"},n_batch,n_annealing,box_length,dopt,n_opt,1);

	// Quit if the MSE goes below 1%
	bmla.set_mse_quit(1.0);

	// Turn on L2 reg
	// bmla.set_l2_reg(0.01);

	// Iterate over the things
	int n_files = 20;
	for (int i=1; i<=n_files; i++) {
		std::cout << "--- File: " << i << " / " << n_files << " ---" << std::endl; 

		// Read the existing as a guess
		// bmla.read("../bimol_annihilation/lattice_v"+pad_str(i,2)+"/lattice/init.txt");		

		// Solve
		bmla.solve("../bimol_annihilation/lattice_v"+pad_str(i,2)+"/lattice/0000.txt",true);
		
		// Write
		bmla.write("../bimol_annihilation/lattice_v"+pad_str(i,2)+"/lattice/init.txt");
	};

	return 0;
};