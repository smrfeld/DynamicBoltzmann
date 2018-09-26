#ifndef MEMORY_H
#define MEMORY_H
#include <memory>
#endif

/************************************
* Namespace for bmla
************************************/

namespace bmla {

	// Forwards and useful typedefs for ixn param
	class IxnParam;

	// Shared ptrs
	typedef std::shared_ptr<IxnParam> Iptr;
};