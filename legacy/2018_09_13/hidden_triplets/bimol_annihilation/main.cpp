#include <iostream>
#include <math.h>

// Get the header
#include "lattice_gillespie.hpp"

using namespace LatticeGillespie;

// Function to get the mean of a vector
double get_mean(std::vector<int> v) {
	double ctr = 0.0;
	for (auto i: v) { ctr += i; };
	return ctr / v.size();
};

// Function to get the std dev of a vector
double get_std(std::vector<int> v) {
	double mean = get_mean(v);
	double ctr = 0.0;
	for (auto i: v) { ctr += pow(i - mean,2); };
	return sqrt(ctr / (v.size()-1));
};

int main() {

	// Seed random no
	srand (time(NULL));
	// srand (2);

	// Box length
	int box_length = 1000;

	// Timestep
	double dt = 0.01;

	// Number of steps to run
	int n_steps = 101;

	// Two initial parameters
	// 100 particles
	// 50 NN
	// Few (1) triplet
	double h_few = -6.85683;
	double j_few = 10.9643; 
	double k_few = -7.96114;
	// Many (40) triplets
	double h_many = -2.99927;
	double j_many = 0.0612041; 
	double k_many = 2.77309;

	// Dict versions
	std::map<std::string,double> h_dict_few,h_dict_many;
	std::map<std::string,std::map<std::string,double>> j_dict_few,j_dict_many;		
	std::map<std::string,std::map<std::string,std::map<std::string,double>>> k_dict_few,k_dict_many;
	h_dict_few["A"] = h_few;
	j_dict_few["A"]["A"] = j_few;
	k_dict_few["A"]["A"]["A"] = k_few;
	h_dict_many["A"] = h_many;
	j_dict_many["A"]["A"] = j_many;
	k_dict_many["A"]["A"]["A"] = k_many;

	// no sampling steps
	int n_sampling = 1000;

	// Indexes to write for mode1, mode2 (inclusive)
	std::map<int,int> mode_min, mode_max;
	mode_min[1] = 1;
	mode_max[1] = 100;
	mode_min[2] = 101;
	mode_max[2] = 200;

	// Two modes
	for (int imode=1; imode<=2; imode++)
	{
		/****************************************
		Make a simulation!
		****************************************/

		for (int write_version=mode_min[imode]; write_version<=mode_max[imode]; write_version++)
		{
			std::cout << "--- Simulating version: " << write_version << " ---" << std::endl;

			// Make a simulation
			Simulation sim(dt,box_length,1); // 1 = 1D

			/********************
			Add species
			********************/

			sim.add_species("A");

			/********************
			Add bimol rxns
			********************/

			// M1
			sim.add_bi_rxn("ann", 0.01, "A","A");

			/********************
			Populate the lattice
			********************/
			
			if (imode==1) {
				sim.populate_lattice(h_dict_few, j_dict_few, k_dict_few, n_sampling);
			} else {
				sim.populate_lattice(h_dict_many, j_dict_many, k_dict_many, n_sampling);
			};

			/********************
			Run
			********************/

			// Run
			bool verbose = false;
			bool write_counts = true;
			bool write_nns = true;
			bool write_latt = true;
			int write_step = 1;
			sim.run(n_steps,verbose,write_counts,write_nns,write_latt,write_step,write_version);

		};
	};

	return 0;
}