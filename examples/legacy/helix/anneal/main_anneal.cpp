#include "dynamicboltz/dynamic_boltzmann.hpp"
#include "dynamicboltz/general.hpp"
#include <iostream>

using namespace DynamicBoltzmann;

int main() {	

	/********************
	Dimensions
	********************/

	std::vector<Dim> dims;
	dims.push_back(Dim("hA",DimType::H,"A",{"hA","hB"},-3.5,0.0,40,-0.8683));
	dims.push_back(Dim("hB",DimType::H,"B",{"hA","hB"},-3.5,0.0,40,-1.3854));
	dims.push_back(Dim("hC",DimType::H,"C",{"hC"},-8.0,0.0,40,-7.0809));

	/********************
	Optimization params
	********************/

	double t_max=1.0;
	int n_t=100;
	int box_length = 1000;

	/********************
	Filenames
	********************/

	std::vector<std::string> fnames;
	fnames.push_back("../stoch_sim/lattice_v001/lattice/");

	/********************
	Options
	********************/

	OptionsSolve options = OptionsSolve();

	// Write only the last basis func
	options.write_bf_only_final = true;

	// Writing
	options.write = true;
	options.dir_write = "data_learned/";

	// Always use the same latt
	options.use_same_lattice_in_batch = true;

	// Exp decay
	options.exp_decay = true;
	for (int i=0; i<200; i++) {
		options.exp_decay_t0_values.push_back(i/20.0);
		options.exp_decay_lambda_values.push_back(0.1);
	};

	/********************
	Solve, annealing
	********************/
 
	// Opt params 
	int batch_size = 10;
	int n_cd_steps = 1;
	double dopt = 1.0;
	int n_opt = 200;

	// Opt problem
	OptProblem opt = OptProblem(dims,{"A","B","C"},{},t_max,n_t,box_length,1);

	// Read in the old basis funcs
	/*
	opt.read_basis_func("F_hA","data_learned/F/F_hA_"+pad_str(i_anneal*n_opt,4)+".txt");
	opt.read_basis_func("F_hB","data_learned/F/F_hB_"+pad_str(i_anneal*n_opt,4)+".txt");
	opt.read_basis_func("F_hC","data_learned/F/F_hC_"+pad_str(i_anneal*n_opt,4)+".txt");
	*/

	// Solve
	std::cout << "Solving..." << std::endl;
	opt.solve(fnames,n_opt,batch_size,n_cd_steps,dopt,options);
	std::cout << "fin." << std::endl;


	return 0;
}