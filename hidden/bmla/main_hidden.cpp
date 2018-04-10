#include <iostream>
#include "bmla.hpp"
#include "general.hpp"
#include <fstream>
#include <vector>

using namespace DynamicBoltzmann;

int main() {

	int box_length = 1000;
	int n_cd_steps = 10;
	double dopt = 0.0002;
	int n_opt = 100;
	int n_batch = 50;

	// Many
	// Few
	// h = -2.0, w = -0.04

	std::vector<Dim> dims;
	dims.push_back(Dim("h",DimType::H,"A",-3.908));
	dims.push_back(Dim("W",DimType::W,"A",1.0078));

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

	// Files to read
	std::vector<std::string> fnames;
	for (int i=1; i<=100; i++) {
		fnames.push_back("../bimol_annihilation/lattice_v"+pad_str(i,3)+"/lattice/0000.txt");
	};

	// Solve
	bmla.solve(fnames,true);
	
	return 0;
};