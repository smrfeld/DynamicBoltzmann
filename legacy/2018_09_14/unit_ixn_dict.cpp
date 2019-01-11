#include "../../include/dynamicboltz_bits/unit_ixn_dict.hpp"

// Other headers
#include "../../include/dynamicboltz_bits/species.hpp"
#include "../../include/dynamicboltz_bits/ixn_param.hpp"

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	BiasDict
	****************************************/

	// Constructor
	BiasDict::BiasDict() {};
	BiasDict::BiasDict(const BiasDict& other) {
		_copy(other);
	};
	BiasDict::BiasDict(BiasDict&& other) {
		_move(other);
	};
    BiasDict& BiasDict::operator=(const BiasDict& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    BiasDict& BiasDict::operator=(BiasDict&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };
	BiasDict::~BiasDict()
	{
		_clean_up();
	};
	void BiasDict::_clean_up() {};
	void BiasDict::_copy(const BiasDict& other) {
		_dict = other._dict;
	};
	void BiasDict::_move(BiasDict& other) {
		_dict = other._dict;
		other._dict.clear();
	};

	// Add to the dict
	void BiasDict::add_ixn(Sptr sp, Iptr ixn) {
		_dict[sp].push_back(ixn);
	};

	// Get from the dict
	double BiasDict::get_ixn(Sptr sp, int timepoint) const {
		auto it1 = _dict.find(sp);
		if (it1 == _dict.end()) {
			return 0.0;
		};

		double r=0.0;
		for (auto ptr: it1->second) {
			r += ptr->get_val_at_timepoint(timepoint);
		};
		return r;
	};
};