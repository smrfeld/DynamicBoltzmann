// string
#ifndef STRING_h
#define STRING_h
#include <string>
#endif

// utility for pair
#ifndef PAIR_h
#define PAIR_h
#include <utility>
#endif

// map
#ifndef MAP_h
#define MAP_h
#include <map>
#endif

/************************************
* Namespace for Gillespie3D
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Species
	****************************************/
	
	struct Species {
		std::string name;
		// Counts
		std::map<Species*,int> nn_count;
		int count;
		// Coupling strengths
		double h;
		std::map<Species*,double> j;
		Species(std::string nameIn);
	};
	// Comparator
	bool operator <(const Species& a, const Species& b);
};

