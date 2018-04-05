#include <iostream>
#include "ixn_param.hpp" // also includes species header

/************************************
* Namespace for Gillespie3D
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Species
	****************************************/
	
	Species::Species(std::string name) {
		_name = name;
		_count = 0;

		// Ptrs
		_h_ptr = nullptr;
		_w_ptr = nullptr;
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
		_nn_count.clear();
		_count = 0;
		_h_ptr = nullptr;
		_j_ptr.clear();
		_w_ptr = nullptr;
	};
	void Species::_copy(const Species& other) {
		_name = other._name;
		_nn_count = other._nn_count;
		_count = other._count;
		_h_ptr = other._h_ptr;
		_j_ptr = other._j_ptr;
		_w_ptr = other._w_ptr;
	};
	/********************
	Set h, j ptr
	********************/

	void Species::set_h_ptr(IxnParam *h_ptr) {
		_h_ptr = h_ptr;
	};
	void Species::add_j_ptr(Species* sp, IxnParam *j_ptr) {
		_j_ptr[sp] = j_ptr;
		// Also add entry in nn count
		_nn_count[sp] = 0;
	};
	void Species::set_w_ptr(IxnParam *w_ptr) {
		_w_ptr = w_ptr;
	};

	/********************
	Setters/getters
	********************/

	double Species::h() const {
		return _h_ptr->get();
	};
	double Species::j(Species *other) const {
		if (_j_ptr.at(other)) {
			return _j_ptr.at(other)->get();
		} else {
			return 0.0;
		};
	};
	double Species::w() const {
		return _w_ptr->get();
	};
	int Species::count() const {
		return _count;
	};
	int Species::nn_count(Species *other) const {
		return _nn_count.at(other);
	};

	std::string Species::name() const { return _name; };

	/********************
	Increment counts
	********************/

	void Species::count_plus() { _count++; };
	void Species::count_minus() { _count--; };
	void Species::nn_count_plus(Species* other) { _nn_count[other]++; };
	void Species::nn_count_minus(Species* other) { _nn_count[other]--; };

	/********************
	Reset counts
	********************/

	void Species::reset_counts() {
		_count = 0;
		for (auto it=_nn_count.begin(); it!=_nn_count.end(); it++) {
			it->second = 0;
		};
	};

	/********************
	Comparator
	********************/

	bool operator <(const Species& a, const Species& b) {
		return a.name() < b.name();
	};

};