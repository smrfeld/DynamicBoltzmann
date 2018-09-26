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
	struct Sptr2;
	struct Sptr3;
	
	// Shared ptrs
	typedef std::shared_ptr<Species> Sptr;
};