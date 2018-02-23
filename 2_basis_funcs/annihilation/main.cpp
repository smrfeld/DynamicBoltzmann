#include "gillespie3d.hpp"
#include <iostream>

using namespace Gillespie3D;

int main() {

	// Seed random no
	srand (2);

	// Initial h, j
	double h0 = 0.5;
	double j0 = 0.1;
	std::map<std::string,double> hdict;
	hdict["A"] = h0;
	std::map<std::string,std::map<std::string,double>> jdict;
	jdict["A"]["A"] = j0;

	// No annealing steps
	int n_annealing = 10000;

	// Alterative: initial counts
	// std::map<std::string,int> counts0;
	// counts0["A"] = 20;

	// Params
	int box_length = 10;
	double dt = 0.01;
	int n_steps = 101;
	bool verbose = false;
	bool write_counts = true;
	bool write_nns = true;
	bool write_latt = true;
	int write_step = 1;

	// Number of samples to run
	int n_samples = 100;

	for (int write_version_no=0; write_version_no<n_samples; write_version_no++)
	{
		std::cout << write_version_no << " / " << n_samples << std::endl;

		Simulation sim(dt,box_length);

		sim.add_species("A");

		sim.add_uni_rxn("rxn", 1.0, "A");

		sim.populate_lattice(hdict,jdict,n_annealing);

		sim.run(n_steps,verbose,write_counts,write_nns,write_latt,write_step,write_version_no);
	};

	return 0;
};