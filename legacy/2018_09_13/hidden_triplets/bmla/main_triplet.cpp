#include <iostream>
#include "bmla.hpp"
#include "general.hpp"
#include <fstream>
#include <vector>

using namespace DynamicBoltzmann;

int main() {

	int box_length = 1000;
	int n_cd_steps = 10; // 100
	double dopt = 0.002;
	int n_opt = 1000;
	int n_batch = 10; // 100

	std::vector<Dim> dims;
	dims.push_back(Dim("h",DimType::H,"A",-2.99927));
	dims.push_back(Dim("J",DimType::J,"A","A",0.0612041));
	dims.push_back(Dim("K",DimType::K,"A","A","A",2.77309));

	BMLA bmla(dims,{"A"},n_batch,box_length,dopt,n_opt,1);

	// Write the solution traj
	bmla.set_write_soln_traj("triplet_soln.txt");

	// Set the number of CD steps
	bmla.set_n_cd_steps(n_cd_steps);

	// Quit if the MSE goes below 1%
	bmla.set_mse_quit(0.1);

	// Turn on L2 reg
	// bmla.set_l2_reg(0.01);

	// Files to read
	std::vector<std::string> fnames;
	for (int i=101; i<=200; i++) {
		fnames.push_back("../bimol_annihilation/lattice_v"+pad_str(i,3)+"/lattice/0000.txt");
	};

	// Solve
	bmla.solve(fnames,true);
	
	return 0;
};