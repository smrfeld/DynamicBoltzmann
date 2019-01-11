#include "dynamicboltz/dynamic_boltzmann.hpp"
#include "dynamicboltz/general.hpp"
#include <iostream>

using namespace DynamicBoltzmann;

int main() {	

	/********************
	Dimensions
	********************/

	std::vector<Dim> dims;
	dims.push_back(Dim("hA",DimType::H,"A",{"hA","hB"},-3.5,0.0,80,-0.8683));
	dims.push_back(Dim("hB",DimType::H,"B",{"hA","hB"},-3.5,0.0,80,-1.3854));
	dims.push_back(Dim("hC",DimType::H,"C",{"hC"},-8.0,0.0,80,-7.0809));

	/********************
	Optimization object
	********************/

	double t_max = 10.0;
	int n_t = 1000;
	int box_length = 1000;

	OptProblem opt = OptProblem(dims,t_max,n_t,box_length,1);

	/********************
	Read in bf
	********************/

	opt.read_basis_func("F_hA","data_learned_split/F/F_hA_0400.txt");
	opt.read_basis_func("F_hB","data_learned_split/F/F_hB_0400.txt");
	opt.read_basis_func("F_hC","data_learned_split/F/F_hC_0400.txt");

	/********************
	Integrate
	********************/

	opt.solve_ixn_param_traj();

	/********************
	Write
	********************/

	opt.write_ixn_params("test_split/", 0);

	return 0;
}