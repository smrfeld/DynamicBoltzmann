#include "species.hpp"

/************************************
* Namespace for Gillespie3D
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Species
	****************************************/
	
	Species::Species(std::string nameIn) {
		this->name = nameIn;
		this->count = 0;
	};

	// Comparator
	bool operator <(const Species& a, const Species& b) {
		return a.name < b.name;
	};
};