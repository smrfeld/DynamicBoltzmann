#include <iostream>
#include "bmla.hpp"
#include "general.hpp"
#include <fstream>
#include <vector>

using namespace DynamicBoltzmann;

int main() {

	// MODE
	// True = Solve for many triplets
	// False = Solve for few triplets
	bool solve_many = true;

	int box_length = 1000;
	int n_cd_steps = 10;
	double dopt = 0.002;
	int n_opt = 10000;
	int n_batch = 10;

	// Initial h,W
	double hinit,winit;
	if (solve_many) {
		// Many
		hinit = 0.475533;
		winit = -3.60267;
	} else {
		// Few
		hinit = 0.860922;
		winit = -5.26262;
	};

	// Dimensions
	std::vector<Dim> dims;
	dims.push_back(Dim("h",DimType::H,"A",hinit));
	dims.push_back(Dim("W",DimType::W,"A",winit));

	// BMLA solver
	BMLA bmla(dims,{"A"},n_batch,box_length,dopt,n_opt,1);

	// Write the solution traj
	if (solve_many) {
		// Many
		bmla.set_write_soln_traj("hidden_soln_many_triplets.txt");
	} else {
		// Few
		bmla.set_write_soln_traj("hidden_soln_few_triplets.txt");
	};

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
	int i_start, i_end;
	if (solve_many) {
		// Many
		i_start = 101;
		i_end = 200;
	} else {
		// Few
		i_start = 1;
		i_end = 100;
	};
	std::vector<std::string> fnames;
	for (int i=i_start; i<=i_end; i++) {
		fnames.push_back("../bimol_annihilation/lattice_v"+pad_str(i,3)+"/lattice/0000.txt");
	};

	// Solve
	bmla.solve(fnames,true);
	
	return 0;
};