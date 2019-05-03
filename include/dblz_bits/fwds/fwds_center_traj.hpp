#ifndef MEMORY_H
#define MEMORY_H
#include <memory>
#endif

/************************************
* Namespace for bmla
************************************/

#ifndef FWDS_CENTER_TRAJ_H
#define FWDS_CENTER_TRAJ_H

namespace dblz {

	// Forwards and useful typedefs for ixn param
	class CenterTraj;

	// Shared ptrs
	typedef std::shared_ptr<CenterTraj> CTptr;
};

#endif
