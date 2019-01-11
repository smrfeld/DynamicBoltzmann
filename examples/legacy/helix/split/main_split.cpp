#include "dynamicboltz/dynamic_boltzmann.hpp"
#include "dynamicboltz/general.hpp"
#include <iostream>

using namespace DynamicBoltzmann;

int main() {	

	/********************
	Dimensions
	********************/

	std::vector<Dim> dims;
	dims.push_back(Dim("hA",DimType::H,"A",{"hA","hB"},-3.5,0.0,80,0.0));
	dims.push_back(Dim("hB",DimType::H,"B",{"hA","hB"},-3.5,0.0,80,0.0));
	dims.push_back(Dim("hC",DimType::H,"C",{"hC"},-8.0,0.0,80,0.0));

	/********************
	Optimization object
	********************/

	double t_max = 0.25;
	int n_t = 25;
	int box_length = 1000;

	OptProblem opt = OptProblem(dims,t_max,n_t,box_length,1);

	/********************
	Filenames
	********************/

	std::vector<std::pair<std::string,int>> fnames;

	/********************
	Options
	********************/

	OptionsSolve options = OptionsSolve();

	// Write only the last basis func
	options.write_bf_only_final = true;

	// Writing
	options.write = true;
	options.dir_write = "data_learned_split/";

	// Fnames to use in every batch
	for (int i=1; i<=8; i++) {
		options.fnames_used_in_every_batch.push_back(std::make_pair("stoch_sim_"+pad_str(i,1)+"/lattice_v001/lattice/",i));
	};

	/********************
	Solve
	********************/

	// Opt params 
	int batch_size = 8;
	int n_cd_steps = 25;
	double dopt = 0.01;
	int n_opt = 400;

	// Solve
	std::cout << "Solving..." << std::endl;
	opt.solve_varying_ic(fnames,n_opt,batch_size,n_cd_steps,dopt,options);
	std::cout << "fin." << std::endl;

	return 0;
}