#include "../include/dblz_bits/fname.hpp"

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

	const std::vector<FName>& FNameColl::get_fnames_all() const {
		return _fnames;
	};
	const FName& FNameColl::get_fname(int idx) const {
		if (idx >= _fnames.size()) {
			std::cerr << ">>> Error: FNameColl::get_fname <<< idx: " << idx << " out of range: " << _fnames.size() << std::endl;
			exit(EXIT_FAILURE);
		};
		return _fnames[idx];
	};
    std::vector<FName> FNameColl::get_fnames_batch(std::vector<int> idxs) const {
        std::vector<FName> fnames;
        for (auto const &idx: idxs) {
            fnames.push_back(get_fname(idx));
        };
        return fnames;
    };
    
	/********************
	Add fname
	********************/

	void FNameColl::add_fname(FName fname) {
		_fnames.push_back(fname);
		_idxs.push_back(_fnames.size()-1);
	};

	/********************
	Get random subset
	********************/

	std::vector<int> FNameColl::get_random_subset(int size) {
		if (_fnames.size() < size) {
			std::cerr << ">>> Error: FNameColl::get_random_subset <<< size of subset is equal to size of filename collection" << std::endl;
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
    
    std::vector<FName> FNameColl::get_random_subset_fnames(int size) {
        std::vector<FName> fnames;
        std::vector<int> idxs = get_random_subset(size);
        for (auto const &idx: idxs) {
            fnames.push_back(_fnames.at(idx));
        };
        return fnames;
    };
};
