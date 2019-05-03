#ifndef MEMORY_H
#define MEMORY_H
#include <memory>
#endif

/************************************
* Namespace for bmla
************************************/

#ifndef FWDS_CENTER_H
#define FWDS_CENTER_H

namespace dblz {

	// Forwards and useful typedefs for ixn param
	class Center;

	// Shared ptrs
	typedef std::shared_ptr<Center> Cptr;
};

#endif
