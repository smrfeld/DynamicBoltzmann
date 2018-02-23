#include "dynamic_boltzmann.hpp"
#include "general.hpp"

using namespace DynamicBoltzmann;

int main() {

	/********************
	Test optimization problem
	********************/

	// h,j
	double h_min = -0.3;
	double h_max = 0.55;
	double h_init = 0.5;
	int n_h = 21;
	double j_min = -0.7;
	double j_max = 0.15;
	double j_init = 0.1;
	int n_j = 21;

	// times
	double t_max=1.0;
	int n_t = 101;

	// Opt params
	int batch_size = 20;
	int n_annealing = 50;
	int box_length = 3;
	double dopt = 0.1;
	int n_opt = 100;
	
	// Init
	OptProblem opt(h_min, h_max, n_h, j_min, j_max, n_j, t_max, n_t, h_init, j_init, batch_size, n_annealing, box_length, dopt, n_opt);

	// Add filenames
	for (int i=0; i<100; i++) {
		opt.add_fname("annihilation/lattice_v" + pad_str(i,2) + "/lattice/");
	};

	// Add species
	opt.add_species("A");

	// Solve
	opt.solve(false);

	return 0;
};