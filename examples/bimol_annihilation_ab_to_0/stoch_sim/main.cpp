#include <iostream>
#include <math.h>

#include <lattGillespie>

using namespace latg;
using namespace std;

int main() {

	// Seed random no
	srand (time(NULL));
	// srand (2);

	// h dict
	map<string,double> h_dict;
	h_dict["A"] = 1.0;
	h_dict["B"] = 1.0;
	map<string,map<string,double>> j_dict;
	j_dict["A"]["A"] = 0.0;
	j_dict["A"]["B"] = 0.0;
	j_dict["B"]["A"] = 0.0;
	j_dict["B"]["B"] = 0.0;

	// Box length
	int box_length = 10;

	// Timestep
	double dt = 0.01;

	// Number of steps to run
	int n_steps = 101;

	// Run
	bool verbose = true;
	bool write_counts = true;
	bool write_nns = false;
	bool write_latt = true;
	int write_step = 1;
	for (int i_version=0; i_version<10; i_version++) {

		// Make a simulation
		Simulation sim(dt,box_length);

		// Add species
		sim.add_species("A");
		sim.add_species("B");

		// Add unimol rxns
		
		// P1
		sim.add_bi_rxn("rxn", 0.1, "A", "B");

		// Populate
		sim.populate_lattice(h_dict,j_dict,1000);

		// Simulate
		sim.run(n_steps,verbose,write_counts,write_nns,write_latt,write_step,i_version);
	};

	return 0;
}