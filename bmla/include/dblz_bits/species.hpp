#include <string>

#include "fwds/fwds_species.hpp"

/************************************
* Namespace for bmla
************************************/

namespace dblz {

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

