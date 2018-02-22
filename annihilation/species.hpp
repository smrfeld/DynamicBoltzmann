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

// Other Gillespie3D

#ifndef REACTIONS_h
#define REACTIONS_h
#include "reactions.hpp"
#endif

/************************************
* Namespace for Gillespie3D
************************************/

namespace Gillespie3D {

	/****************************************
	Necessary declarations
	****************************************/

	struct Mol;

	/****************************************
	Species
	****************************************/
	
	struct Species {
		std::string name;
		std::map<Species*,std::vector<BiReaction*>> bi_rxns;
		std::vector<UniReaction*> uni_rxns;
		bool conserved;
		int count;
		Species(std::string nameIn, bool conservedIn);

		/********************
		Add a reaction, if appropriate
		********************/

		void add_rxn(BiReaction* rxn);
		void add_rxn(UniReaction* rxn);

		/********************
		Check if any reactions are possible; if so, return a random one
		********************/

		std::pair<bool,BiReaction*> check_bi_rxns_mol(Mol *other);

	};	
	// Comparator
	bool operator <(const Species& a, const Species& b);

	/****************************************
	Mol
	****************************************/
	
	struct Mol {
		Species *sp;
		Mol(Species *spIn);

		/********************
		Check if any reactions are possible; if so, return a random one
		********************/

		std::pair<bool,BiReaction*> check_bi_rxns_mol(Mol *other);
	};
	// Comparator
	bool operator <(const Mol& a, const Mol& b);
};