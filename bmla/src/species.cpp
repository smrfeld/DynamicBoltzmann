#include <iostream>
#include "ixn_param.hpp" // also includes species header

/************************************
* Namespace for Gillespie3D
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Doublets, Triplets of Species ptrs
	****************************************/

	Species2::Species2(Species* sp1, Species* sp2) {
		if (sp1>=sp2) {
			this->sp1 = sp1;
			this->sp2 = sp2;
		} else {
			this->sp1 = sp2;
			this->sp2 = sp1;
		};
	};
	bool operator <(const Species2& a, const Species2& b) {
    	return std::tie(a.sp1, a.sp2) < std::tie(b.sp1, b.sp2);
	};

	Species3::Species3(Species* sp1, Species* sp2, Species *sp3) {
		if (sp1>=sp2 && sp2 >= sp3) {
			this->sp1 = sp1;
			this->sp2 = sp2;
			this->sp3 = sp3;
		} else if (sp1>=sp3 && sp3 >= sp2) {
			this->sp1 = sp1;
			this->sp2 = sp3;
			this->sp3 = sp2;
		} else if (sp2>=sp1 && sp1 >= sp3) {
			this->sp1 = sp2;
			this->sp2 = sp1;
			this->sp3 = sp3;
		} else if (sp2>=sp3 && sp3 >= sp1) {
			this->sp1 = sp2;
			this->sp2 = sp3;
			this->sp3 = sp1;
		} else if (sp3>=sp1 && sp1 >= sp2) {
			this->sp1 = sp3;
			this->sp2 = sp1;
			this->sp3 = sp2;
		} else if (sp3>=sp2 && sp2 >= sp1) {
			this->sp1 = sp3;
			this->sp2 = sp2;
			this->sp3 = sp1;
		};
	};
	bool operator <(const Species3& a, const Species3& b) {
    	return std::tie(a.sp1, a.sp2, a.sp3) < std::tie(b.sp1, b.sp2, b.sp3);
	};

	/****************************************
	Species
	****************************************/

	Species::Species(std::string name) {
		_name = name;

		// Counters
		_count = nullptr;

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
		_count = nullptr;
		_nn_count.clear();
		_triplet_count.clear();
		_quartic_count.clear();
		_h_ptr = nullptr;
		_j_ptr.clear();
		_k_ptr.clear();
		_w_ptr = nullptr;
	};
	void Species::_copy(const Species& other) {
		_name = other._name;
		_count = other._count;
		_nn_count = other._nn_count;
		_triplet_count = other._triplet_count;
		_quartic_count = other._quartic_count;
		_h_ptr = other._h_ptr;
		_j_ptr = other._j_ptr;
		_k_ptr = other._k_ptr;
		_w_ptr = other._w_ptr;
	};

	/********************
	Set counters
	********************/

	void Species::set_counter(Counter *ctr) {
		_count = ctr;
	};
	void Species::add_nn_counter(Species *sp, Counter *ctr) {
		_nn_count[sp] = ctr;
	};
	void Species::add_triplet_counter(Species *sp1, Species *sp2, Counter *ctr) {
		_triplet_count[Species2(sp1,sp2)] = ctr;
	};
	void Species::add_quartic_counter(Species *sp1, Species *sp2, Species *sp3, Counter *ctr) {
		_quartic_count[Species3(sp1,sp2,sp3)] = ctr;
	};

	/********************
	Set h, j ptr
	********************/

	void Species::set_h_ptr(IxnParam *h_ptr) {
		_h_ptr = h_ptr;
	};
	void Species::add_j_ptr(Species* sp, IxnParam *j_ptr) {
		_j_ptr[sp] = j_ptr;
	};
	void Species::add_k_ptr(Species* sp1, Species* sp2, IxnParam *k_ptr) {
		_k_ptr[Species2(sp1,sp2)] = k_ptr;
	};
	void Species::set_w_ptr(IxnParam *w_ptr) {
		_w_ptr = w_ptr;
	};

	/********************
	Get ixn params
	********************/

	double Species::h() const {
		return _h_ptr->get();
	};
	double Species::j(Species *other) const {
		auto it = _j_ptr.find(other);
		if (it != _j_ptr.end()) {
			return it->second->get();
		} else {
			return 0.0;
		};
	};
	double Species::w() const {
		return _w_ptr->get();
	};
	double Species::k(Species* other1, Species *other2) const {
		auto it = _k_ptr.find(Species2(other1,other2));
		if (it != _k_ptr.end()) {
			return it->second->get();
		} else {
			return 0.0;
		};
	};

	/********************
	Get counts
	********************/

	double Species::count() const {
		return _count->get_count();
	};
	double Species::nn_count(Species *other) const {
		return _nn_count.at(other)->get_count();
	};
	double Species::triplet_count(Species* other1, Species *other2) const {
		return _triplet_count.at(Species2(other1,other2))->get_count();
	};
	double Species::quartic_count(Species* other1, Species *other2, Species *other3) const {
		return _quartic_count.at(Species3(other1,other2,other3))->get_count();
	};

	/********************
	Name
	********************/

	std::string Species::name() const { return _name; };

	/********************
	Increment counts
	********************/

	void Species::count_increment(double inc) {
		_count->increment(inc);
	};
	void Species::nn_count_increment(Species* other, double inc) {
		_nn_count.at(other)->increment(inc);
	};
	void Species::triplet_count_increment(Species* other1, Species* other2, double inc) {
		_triplet_count.at(Species2(other1,other2))->increment(inc);
	};
	void Species::quartic_count_increment(Species* other1, Species* other2, Species *other3, double inc) {
		_quartic_count.at(Species3(other1,other2,other3))->increment(inc);
	};

	/********************
	Reset counts
	********************/

	void Species::reset_counts() {
		_count->reset_count();
		for (auto it=_nn_count.begin(); it!=_nn_count.end(); it++) {
			it->second->reset_count();
		};
		for (auto it1=_triplet_count.begin(); it1!=_triplet_count.end(); it1++) {
			it1->second->reset_count();
		};
		for (auto it1=_quartic_count.begin(); it1!=_quartic_count.end(); it1++) {
			it1->second->reset_count();
		};
	};

	/********************
	Comparator
	********************/

	bool operator <(const Species& a, const Species& b) {
		return a.name() < b.name();
	};

};