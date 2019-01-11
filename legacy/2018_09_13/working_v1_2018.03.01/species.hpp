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
		// Name
		std::string name;
		
		// Counts
		std::map<Species*,int> nn_count;
		int count;

		// Pointers...
		// Solution array
		double ***_soln_traj_ptr;
		// Current time in the optimization
		int *_t_opt_ptr;
		int _h_index;
		std::map<Species*,int> _j_index;

		// Constructor
		Species(std::string nameIn);

		// Accessor h,j
		double h();
		double j(Species *other);
	};
	// Comparator
	bool operator <(const Species& a, const Species& b);
};

