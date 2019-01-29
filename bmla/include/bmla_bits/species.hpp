#include <string>

#ifndef FWDS_SPECIES_H
#define FWDS_SPECIES_H
#include "fwds/fwds_species.hpp"
#endif

/************************************
* Namespace for bmla
************************************/

namespace bmla {

	/****************************************
	Species
	****************************************/

	class Species {

	private:

		// Name
		std::string _name;

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const Species& other);

	public:

		/********************
		Constructor
		********************/

		Species(std::string name);
		Species(const Species& other);
		Species(Species&& other);
		Species& operator=(const Species& other);
		Species& operator=(Species&& other);
		~Species();

		/********************
		Name
		********************/

		std::string get_name() const;
	};
	// Comparator
	bool operator <(const Species& a, const Species& b);
};

