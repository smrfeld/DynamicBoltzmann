#include <iostream>
#include "bmla.hpp"
#include "general.hpp"
#include <fstream>
#include <vector>

using namespace DynamicBoltzmann;

int main() {

	// MODE
	// True = Solve for many triplets
	// False = Solve for few triplets
	bool solve_many = false;

	int box_length = 1000;

	// Initial h,W
	double hinit,winit,binit;
	if (solve_many) {
		// Many
		hinit = -0.8388;
		winit = 0.00762406;
		binit = -4.99998;
	} else {
		// Few
		hinit = -0.86438;
		winit = -1.8263;
		binit = -4.96677;
	};

	// Dimensions
	std::vector<Dim> dims;
	dims.push_back(Dim("h",DimType::H,"A",hinit));
	dims.push_back(Dim("W",DimType::W,"A",winit));
	dims.push_back(Dim("B",DimType::B,"A",binit));

	// BMLA solver
	BMLA bmla(dims,{"A"},0,box_length,0,0,1);

	// Add hidden units

	// Connect every visible unit to a hidden unit
	// Each hidden unit is connected to 3 visible units, except at the ends where we only connect to 2
	std::vector<int> lattice_idxs;

	// End #1
	lattice_idxs.clear();
	lattice_idxs.push_back(1);
	lattice_idxs.push_back(2);
	bmla.add_hidden_unit(lattice_idxs, "A");
	// End #2
	lattice_idxs.clear();
	lattice_idxs.push_back(box_length-1);
	lattice_idxs.push_back(box_length);
	bmla.add_hidden_unit(lattice_idxs, "A");
	// Everything in the middle
	for (int idx = 2; idx<=box_length-1; idx++) {
		lattice_idxs.clear();
		lattice_idxs.push_back(idx-1);
		lattice_idxs.push_back(idx);
		lattice_idxs.push_back(idx+1);
		bmla.add_hidden_unit(lattice_idxs, "A");
	};

	// Sample
	bmla.sample(1000,1000,true,true,true);
	
	return 0;
};