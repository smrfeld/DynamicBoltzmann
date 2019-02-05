#ifndef MEMORY_H
#define MEMORY_H
#include <memory>
#endif

/************************************
* Namespace for bmla
************************************/

#ifndef FWDS_IXN_PARAM_H
#define FWDS_IXN_PARAM_H

namespace dblz {

	// Forwards and useful typedefs for ixn param
	class IxnParam;

	// Shared ptrs
	typedef std::shared_ptr<IxnParam> Iptr;
};

#endif
