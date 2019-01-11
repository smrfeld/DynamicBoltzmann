#include <iostream>
#include "bmla/bmla.hpp"
#include "bmla/general.hpp"
#include <fstream>
#include <vector>

using namespace DynamicBoltzmann;

int main() {

	// Box length
	int box_length = 1000;

	/********************
	Dimensions
	********************/

	std::vector<Dim> dims;
	dims.push_back(Dim("hA",DimType::H,"A",0.0));
	dims.push_back(Dim("hB",DimType::H,"B",0.0));
	dims.push_back(Dim("hC",DimType::H,"C",0.0));

	/********************
	Solver object
	********************/

	// BMLA solver
	BMLA bmla(dims,box_length,1);

	/********************
	Files to read
	********************/

	std::vector<std::string> fnames;
	fnames.push_back("stoch_sim/lattice_v001/lattice/0000.txt");

	/********************
	Options
	********************/

	OptionsSolveBMLA options = OptionsSolveBMLA();

	// Write the solution traj
	options.write_soln_traj = true;
	options.fname_write_soln_traj = "ic_params.txt";
	options.write_moment_traj = true;
	options.fname_write_moment_traj = "ic_moments.txt";

	// Solve
	int n_cd_steps = 1;
	double dopt = 0.001;
	int n_opt = 2000;
	int batch_size = 10;
	bmla.solve(fnames, n_opt, batch_size, n_cd_steps, dopt, options);
	
	return 0;
};