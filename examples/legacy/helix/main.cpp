#include "dynamicboltz/dynamic_boltzmann.hpp"
#include "dynamicboltz/general.hpp"
#include <iostream>

using namespace DynamicBoltzmann;

int main() {	

	/********************
	Dimensions
	********************/

	std::vector<Dim> dims;
	dims.push_back(Dim("hA",DimType::H,"A",{"hA","hB"},-3.5,0.0,20,-0.8683));
	dims.push_back(Dim("hB",DimType::H,"B",{"hA","hB"},-3.5,0.0,20,-1.3854));
	dims.push_back(Dim("hC",DimType::H,"C",{"hC"},-8.0,0.0,20,-7.0809));

	/********************
	Optimization object
	********************/

	double t_max = 2.0;
	int n_t = 200;
	int box_length = 1000;

	OptProblem opt = OptProblem(dims,t_max,n_t,box_length,1);

	/********************
	Filenames
	********************/

	std::vector<std::string> fnames;
	fnames.push_back("stoch_sim/lattice_v001/lattice/");

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

	// Local decay factor
	//options.local_decay = true;
	//options.local_decay_factor = 0.1;

	// Write var terms
	// options.write_var_terms = true;

	/********************
	Solve
	********************/

	// Opt params 
	int batch_size = 1;
	int n_cd_steps = 10;
	double dopt = 0.1;
	int n_opt = 400;

	// Solve
	std::cout << "Solving..." << std::endl;
	opt.solve(fnames,n_opt,batch_size,n_cd_steps,dopt,options);
	std::cout << "fin." << std::endl;

	return 0;
}