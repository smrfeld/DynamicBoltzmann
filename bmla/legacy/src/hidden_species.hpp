// string
#ifndef STRING_h
#define STRING_h
#include <string>
#endif

/************************************
* Namespace for Gillespie3D
************************************/

namespace DynamicBoltzmann {

	/****************************************
	HiddenSpecies
	****************************************/
	
	class HiddenSpecies {

	private:

		// Name
		std::string _name;

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const HiddenSpecies& other);

	public:

		// Constructor
		HiddenSpecies(std::string name);
		HiddenSpecies(const HiddenSpecies& other);
		HiddenSpecies(HiddenSpecies&& other);
		HiddenSpecies& operator=(const HiddenSpecies& other);
		HiddenSpecies& operator=(HiddenSpecies&& other);
		~HiddenSpecies();

		// Name
		std::string name() const;
	};
	// Comparator
	bool operator <(const HiddenSpecies& a, const HiddenSpecies& b);
};

