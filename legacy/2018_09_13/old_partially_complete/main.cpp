#include "dynamic_boltzmann.hpp"

using namespace DynamicBoltzmann;

int main() {

	/********************
	Test interpolation
	********************/

	BasisFunc bf(0.0,1.0,11);
	bf.sin_func();
	bf.update_derivs();
	// std::cout << l << std::endl;

	bf.write_to_file("l1.txt");
	bf.write_to_file("l2.txt",1000);
	bf.write_deriv_to_file("l3.txt");
	bf.write_deriv_to_file("l4.txt",1000);

	/********************
	Test optimization
	********************/

	OptProblem opt(0.0, 1.0, 11, 1.0, 11, 0.3);
	opt.solve_nu_traj();
	//opt.print_nu_traj();
	opt.solve_var_traj();
	opt.print_var_traj();

	return 0;
};