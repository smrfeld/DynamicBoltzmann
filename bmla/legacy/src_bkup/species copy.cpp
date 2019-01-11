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
		//_count = 0;

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
		/*
		_count = 0;
		_nn_count.clear();
		_triplet_count.clear();
		_quartic_count.clear();
		*/
		_counters.clear();
		_h_ptr = nullptr;
		_j_ptr.clear();
		_w_ptr = nullptr;
		_k_ptr.clear();
	};
	void Species::_copy(const Species& other) {
		_name = other._name;
		/*
		_count = other._count;
		_nn_count = other._nn_count;
		_triplet_count = other._triplet_count;
		_quartic_count = other._quartic_count;
		*/
		_counters = other._counters;
		_h_ptr = other._h_ptr;
		_j_ptr = other._j_ptr;
		_w_ptr = other._w_ptr;
		_k_ptr = other._k_ptr;
	};

	/********************
	Add a counter
	********************/

	void Species::add_counter(Counter *counter) {
		_counters.push_back(counter);
	};

	/********************
	Initialize counts by informing this species of the existence of all others
	********************/

	/*
	void Species::init_count_nns(std::list<Species>& sp_list) {
		for (auto sp_it=sp_list.begin(); sp_it != sp_list.end(); sp_it++) {
			init_count_nns(&(*sp_it));
		};
	};
	void init_count_nns(Species* sp);
	void init_count_triplets(std::list<Species>& sp_list);
	void init_count_triplets(Species* sp);
	void init_count_quartics(std::list<Species>& sp_list);
	void init_count_quartics(Species* sp);

	void Species::init_counts(std::list<Species>& sp_list) {
		// Make all the entries
		for (auto sp_it=sp_list.begin(); sp_it != sp_list.end(); sp_it++) {
			_nn_count[&(*sp_it)] = 0.0;
			for (auto sp_it2=sp_list.begin(); sp_it2 != sp_list.end(); sp_it2++) {
				_triplet_count[&(*sp_it)][&(*sp_it2)] = 0.0;
			};	
		};
	};
	*/

	/********************
	Set h, j ptr
	********************/

	void Species::set_h_ptr(IxnParam *h_ptr) {
		_h_ptr = h_ptr;
	};
	void Species::add_j_ptr(Species* sp, IxnParam *j_ptr) {
		_j_ptr[sp] = j_ptr;
	};
	void Species::set_w_ptr(IxnParam *w_ptr) {
		_w_ptr = w_ptr;
	};
	void Species::add_k_ptr(Species* sp1, Species* sp2, IxnParam *k_ptr) {
		_k_ptr[sp1][sp2] = k_ptr;
		_k_ptr[sp2][sp1] = k_ptr;
	};

	/********************
	Setters/getters
	********************/

	double Species::h() const {
		return _h_ptr->get();
	};
	double Species::j(Species *other) const {
		if (_j_ptr.at(other)) { // Check against nullptr
			return _j_ptr.at(other)->get();
		} else {
			return 0.0;
		};
	};
	double Species::w() const {
		return _w_ptr->get();
	};
	double Species::k(Species* other1, Species *other2) const {
		if (_k_ptr.at(other1).at(other2)) { // Check against nullptr
			return _k_ptr.at(other1).at(other2)->get();
		} else {
			return 0.0;
		};
	};
	/*
	double Species::count() const {
		return _count;
	};
	double Species::nn_count(Species *other) const {
		return _nn_count.at(other);
	};
	double Species::triplet_count(Species* other1, Species *other2) const {
		return _triplet_count.at(other1).at(other2);
	};
	*/

	std::string Species::name() const { return _name; };

	/********************
	Increment counts
	********************/

	void Species::counts_increment(double inc) {
		for (auto c: _counters) {
			c->increment(inc);
		};
	};
	void Species::counts_plus() {
		for (auto c: _counters) {
			c->plus();
		};
	};
	void Species::counts_minus() {
		for (auto c: _counters) {
			c->minus();
		};
	};

	/*
	void Species::count_increment(double inc) { 
		_count += inc; };
	void Species::nn_count_increment(Species* other, double inc) { 
		_nn_count[other] += inc; 
	};
	void Species::triplet_count_increment(Species* other1, Species* other2, double inc) {
		_triplet_count[other1][other2] += inc; 
		if (other1 != other2) {
			_triplet_count[other2][other1] += inc; 
		};
	};
	void Species::quartic_count_increment(Species* other1, Species* other2, Species *other3, double inc) {
		_quartic_count[other1][other2][other3] += inc;
		if (other1 != other2) {
			_quartic_count[other2][other1][other3] += inc; 
		};
		if (other2 != other3) {
			_quartic_count[other1][other3][other2] += inc; 
		};
		if (other2 != other3 && other2 != other1) {
			_quartic_count[other2][other3][other1] += inc; 
		};
		if (other1 != other3) {
			_quartic_count[other3][other2][other1] += inc; 
		};
		if (other2 != other3 && other2 != other1) {
			_quartic_count[other3][other1][other2] += inc; 
		};
	};
	*/

	/********************
	Reset counts
	********************/

	void Species::reset_counts() {
		_count = 0.0;
		for (auto it=_nn_count.begin(); it!=_nn_count.end(); it++) {
			it->second = 0.0;
		};
		for (auto it1=_triplet_count.begin(); it1!=_triplet_count.end(); it1++) {
			for (auto it2=it1->second.begin(); it2!=it1->second.end(); it2++) {
				it2->second = 0.0;
			};
		};
	};

	/********************
	Comparator
	********************/

	bool operator <(const Species& a, const Species& b) {
		return a.name() < b.name();
	};

};