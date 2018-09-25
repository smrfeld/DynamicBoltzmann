#include <string>

#ifndef FWDS_SPECIES_H
#define FWDS_SPECIES_H
#include "fwds/fwds_species.hpp"
#endif

/************************************
* Namespace for dblz
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






























	/****************************************
	Doublets, Triplets of Species ptrs
	****************************************/

	// Two species
	struct Sptr2 {
		Sptr s1;
		Sptr s2;
		Sptr2(Sptr s1, Sptr s2);
	};
	bool operator <(const Sptr2& a, const Sptr2& b);

	// Three species
	struct Sptr3 {
		Sptr s1;
		Sptr s2;
		Sptr s3;
		Sptr3(Sptr s1, Sptr s2, Sptr s3);
	};
	bool operator <(const Sptr3& a, const Sptr3& b);
};

