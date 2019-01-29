#ifndef MEMORY_H
#define MEMORY_H
#include <memory>
#endif

/************************************
* Namespace for bmla
************************************/

namespace bmla {

	// Forwards and useful typedefs for species
	class Species;
	
	// Shared ptrs
	typedef std::shared_ptr<Species> Sptr;
};
