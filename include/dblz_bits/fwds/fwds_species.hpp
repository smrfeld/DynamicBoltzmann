#ifndef MEMORY_H
#define MEMORY_H
#include <memory>
#endif

/************************************
* Namespace for bmla
************************************/

#ifndef FWDS_SPECIES_H
#define FWDS_SPECIES_H

namespace dblz {

	// Forwards and useful typedefs for species
	class Species;
	
	// Shared ptrs
	typedef std::shared_ptr<Species> Sptr;
};

#endif
