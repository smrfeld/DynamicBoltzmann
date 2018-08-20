#include "species.hpp"

// Other headers
#include "counter.hpp"
#include "ixn_param_traj.hpp"

#include <iostream>

/************************************
* Namespace for dboltz
************************************/

namespace dboltz {

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

	SpeciesVH::SpeciesVH(Species* sp_visible, HiddenSpecies *sp_hidden) {
		this->sp_visible = sp_visible;
		this->sp_hidden = sp_hidden;
	};
	bool operator <(const SpeciesVH& a, const SpeciesVH& b) {
    	return std::tie(a.sp_visible, a.sp_hidden) < std::tie(b.sp_visible, b.sp_hidden);
	};

	/****************************************
	Species
	****************************************/

	Species::Species(std::string name) {
		_name = name;

		// Counters
		_count = nullptr;
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
		_h_ptrs.clear();
		_j_ptrs.clear();
		_k_ptrs.clear();
	};
	void Species::_copy(const Species& other) {
		_name = other._name;
		_count = other._count;
		_nn_count = other._nn_count;
		_triplet_count = other._triplet_count;
		_quartic_count = other._quartic_count;
		_h_ptrs = other._h_ptrs;
		_j_ptrs = other._j_ptrs;
		_k_ptrs = other._k_ptrs;
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

	void Species::add_h_ptr(IxnParamTraj *h_ptr) {
		// Check not already entered
		auto it = std::find(_h_ptrs.begin(),_h_ptrs.end(), h_ptr);
		if (it == _h_ptrs.end()) {
			_h_ptrs.push_back(h_ptr);
		};
	};
	void Species::add_j_ptr(Species* sp, IxnParamTraj *j_ptr) {
		// Check not already entered
		auto it1 = _j_ptrs.find(sp);
		if (it1 == _j_ptrs.end()) {
			_j_ptrs[sp].push_back(j_ptr);
		} else {
			auto it2 = std::find(_j_ptrs[sp].begin(),_j_ptrs[sp].end(), j_ptr);
			if (it2 == _j_ptrs[sp].end()) {
				_j_ptrs[sp].push_back(j_ptr);
			};
		};
	};
	void Species::add_k_ptr(Species* sp1, Species* sp2, IxnParamTraj *k_ptr) {
		Species2 sp = Species2(sp1,sp2);
		// Check not already entered
		auto it1 = _k_ptrs.find(sp);
		if (it1 == _k_ptrs.end()) {
			_k_ptrs[sp].push_back(k_ptr);
		} else {
			auto it2 = std::find(_k_ptrs[sp].begin(),_k_ptrs[sp].end(), k_ptr);
			if (it2 == _k_ptrs[sp].end()) {
				_k_ptrs[sp].push_back(k_ptr);
			};
		};
	};

	/********************
	Get ixn params
	********************/

	double Species::h() const {
		double act=0.0;
		for (auto h: _h_ptrs) {
			act += h->get();
		};
		return act;
	};
	double Species::j(Species *other) const {
		double act=0.0;
		auto it = _j_ptrs.find(other);
		if (it != _j_ptrs.end()) {
			for (auto j: it->second) {
				act += j->get();
			};
		};
		return act;
	};
	double Species::k(Species* other1, Species *other2) const {
		double act=0.0;
		auto it = _k_ptrs.find(Species2(other1,other2));
		if (it != _k_ptrs.end()) {
			for (auto k: it->second) {
				act += k->get();
			};
		};
		return act;
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