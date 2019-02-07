#include "../include/dblz_bits/fname_traj.hpp"

#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>

/************************************
* Namespace for bmla
************************************/

namespace dblz {

	/****************************************
	Filename
	****************************************/

	FName::FName(std::string name, bool binary) {
		this->name = name;
		this->binary = binary;
	};

	/****************************************
	Filename collection
	****************************************/

	/********************
	Get fnames
	********************/

	const std::vector<FNameTraj>& FNameTrajColl::get_fname_trajs_all() const {
		return _fnames;
	};
	const FNameTraj& FNameTrajColl::get_fname_traj(int idx) const {
		if (idx >= _fnames.size()) {
			std::cerr << ">>> Error: FNameTrajColl::get_fname_traj <<< idx: " << idx << " out of range: " << _fnames.size() << std::endl;
			exit(EXIT_FAILURE);
		};
		return _fnames[idx];
	};
    std::vector<FNameTraj> FNameTrajColl::get_fname_trajs_batch(std::vector<int> idxs) const {
        std::vector<FNameTraj> fnames;
        for (auto const &idx: idxs) {
            fnames.push_back(get_fname_traj(idx));
        };
        return fnames;
    };
    
	/********************
	Add fname
	********************/

	void FNameTrajColl::add_fname_traj(FNameTraj fname) {
		_fnames.push_back(fname);
		_idxs.push_back(_fnames.size()-1);
	};

	/********************
	Get random subset
	********************/

	std::vector<int> FNameTrajColl::get_random_subset(int size) {
		if (_fnames.size() < size) {
			std::cerr << ">>> Error: FNameTrajColl::get_random_subset <<< size of subset is equal to size of filename collection" << std::endl;
			exit(EXIT_FAILURE);
		};

		// Shuffle idxs
	    std::random_device rd;
	    std::mt19937 g(rd());
	    std::shuffle(_idxs.begin(), _idxs.end(), g);

	    std::vector<int>::const_iterator first = _idxs.begin();
	    std::vector<int>::const_iterator last = _idxs.begin() + size;
	    return std::vector<int>(first,last);
	};
};
