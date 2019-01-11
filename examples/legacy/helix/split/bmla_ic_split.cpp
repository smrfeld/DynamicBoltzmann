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

	for (int i=1; i<=8; i++)
	{

		std::vector<std::string> fnames;
		fnames.push_back("stoch_sim_"+pad_str(i,1)+"/lattice_v001/lattice/0000.txt");

		/********************
		Options
		********************/

		OptionsSolveBMLA options = OptionsSolveBMLA();

		// Solve
		int n_cd_steps = 1;
		double dopt = 0.001;
		int n_opt = 2000;
		int batch_size = 10;
		bmla.solve(fnames, n_opt, batch_size, n_cd_steps, dopt, options);
		
		// Write
		bmla.write("stoch_sim_"+pad_str(i,1)+"/lattice_v001/init.txt");

	};

	return 0;
};