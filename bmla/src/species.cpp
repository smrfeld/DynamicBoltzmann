#include "../include/bmla_bits/species.hpp"

#include <iostream>

/************************************
* Namespace for bmla
************************************/

namespace bmla {

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
};
