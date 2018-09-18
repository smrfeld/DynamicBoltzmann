#include "../../include/dynamicboltz_bits/species.hpp"

#include <iostream>

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	Species
	****************************************/

	Species::Species(std::string name) {
		_name = name;
	};
	Species::Species(const Species& other) {
		_copy(other);
	};
	Species::Species(Species&& other) {
		_copy(other);
		other._reset();
	};
	Species& Species::operator=(const Species& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	Species& Species::operator=(Species&& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
	};
	Species::~Species() {
		_clean_up();
	};

	void Species::_clean_up() {
		// Nothing...
	};
	void Species::_reset() {
		_name = "";	
	};
	void Species::_copy(const Species& other) {
		_name = other._name;
	};

	/********************
	Name
	********************/

	std::string Species::get_name() const { return _name; };

	/********************
	Comparator
	********************/

	bool operator <(const Species& a, const Species& b) {
		return a.get_name() < b.get_name();
	};




































	/****************************************
	Doublets, Triplets of Species ptrs
	****************************************/

	// Two species
	Sptr2::Sptr2(Sptr s1, Sptr s2) {
		if (s1 >= s2) {
			this->s1 = s1;
			this->s2 = s2;
		} else {
			this->s1 = s2;
			this->s2 = s1;
		};
	};
	bool operator <(const Sptr2& a, const Sptr2& b) {
    	return std::tie(a.s1, a.s2) < std::tie(b.s1, b.s2);
	};

	// Three species
	Sptr3::Sptr3(Sptr s1, Sptr s2, Sptr s3) {
		if (s1>=s2 && s2 >= s3) {
			this->s1 = s1;
			this->s2 = s2;
			this->s3 = s3;
		} else if (s1>=s3 && s3 >= s2) {
			this->s1 = s1;
			this->s2 = s3;
			this->s3 = s2;
		} else if (s2>=s1 && s1 >= s3) {
			this->s1 = s2;
			this->s2 = s1;
			this->s3 = s3;
		} else if (s2>=s3 && s3 >= s1) {
			this->s1 = s2;
			this->s2 = s3;
			this->s3 = s1;
		} else if (s3>=s1 && s1 >= s2) {
			this->s1 = s3;
			this->s2 = s1;
			this->s3 = s2;
		} else if (s3>=s2 && s2 >= s1) {
			this->s1 = s3;
			this->s2 = s2;
			this->s3 = s1;
		};
	};
	bool operator <(const Sptr3& a, const Sptr3& b) {
    	return std::tie(a.s1, a.s2, a.s3) < std::tie(b.s1, b.s2, b.s3);
	};
};