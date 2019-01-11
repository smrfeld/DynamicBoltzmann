#include <iostream>
#include "hidden_species.hpp"

/************************************
* Namespace for Gillespie3D
************************************/

namespace DynamicBoltzmann {

	/****************************************
	HiddenSpecies
	****************************************/

	HiddenSpecies::HiddenSpecies(std::string name) {
		_name = name;
	};
	HiddenSpecies::HiddenSpecies(const HiddenSpecies& other) {
		_copy(other);
	};
	HiddenSpecies::HiddenSpecies(HiddenSpecies&& other) {
		_copy(other);
		other._reset();
	};
	HiddenSpecies& HiddenSpecies::operator=(const HiddenSpecies& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	HiddenSpecies& HiddenSpecies::operator=(HiddenSpecies&& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
	};
	HiddenSpecies::~HiddenSpecies() {
		_clean_up();
	};

	void HiddenSpecies::_clean_up() {
		// Nothing...
	};
	void HiddenSpecies::_reset() {
		_name = "";	
	};
	void HiddenSpecies::_copy(const HiddenSpecies& other) {
		_name = other._name;
	};

	/********************
	Name
	********************/

	std::string HiddenSpecies::name() const { return _name; };

	/********************
	Comparator
	********************/

	bool operator <(const HiddenSpecies& a, const HiddenSpecies& b) {
		return a.name() < b.name();
	};

};