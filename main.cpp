#include "dynamic_boltzmann.hpp"

using namespace DynamicBoltzmann;

int main() {

	/********************
	Test interpolation
	********************/
	/*

	BasisFunc bf(0.0,1.0,11);
	bf.test_fill_sin();
	// std::cout << l << std::endl;

	bf.write_to_file("l1.txt");
	bf.write_to_file("l2.txt",1000);
	bf.write_deriv_to_file("l3.txt");
	bf.write_deriv_to_file("l4.txt",1000);
	*/

	/********************
	Test solving nu, var trajs
	********************/
	/*
	double nu_min=0.0;
	double nu_max=1.0;
	int n_nu=11;
	double t_max=2.0;
	int n_t = 21;
	double nu_init = 0.3;
	int batch_size = 10;
	int n_annealing = 100;
	std::vector<std::string> species_list;
	species_list.push_back("S");
	int box_length = 5;
	double dopt = 0.00001;

	OptProblem opt(nu_min, nu_max, n_nu, t_max, n_t, nu_init, batch_size, n_annealing, species_list, box_length, dopt);
	opt.solve_nu_traj();
	opt.solve_var_traj();
	opt.write_nu_traj("nu_traj.txt");
	opt.write_var_traj("var_traj.txt");
	*/

	/********************
	Test reading lattice
	********************/

	/*
	Lattice l(10);
	l.read_from_file("latt_st_test.txt");
	*/

	/********************
	Test optimization problem
	********************/

	double nu_min=0.0;
	double nu_max=1.0;
	int n_nu=11;
	double t_max=2.0;
	int n_t = 21;
	double nu_init = 0.3;
	int batch_size = 10;
	int n_annealing = 100;
	std::vector<std::string> species_list;
	species_list.push_back("A");
	int box_length = 20;
	double dopt = 0.00001;

	OptProblem opt(nu_min, nu_max, n_nu, t_max, n_t, nu_init, batch_size, n_annealing, species_list, box_length, dopt);

	// Write the initial bfs
	opt.write_bfs("f0.txt");

	// Solve
	opt.solve();

	// Write the final bfs
	opt.write_bfs("f1.txt");

	return 0;
};