#include <iostream>
#include "ixn_param_traj.hpp" // also includes species header

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
		_t_opt_ptr = nullptr;
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
		_triplet_count.clear();
		_count = 0;
		_t_opt_ptr = nullptr;
		_h_ptr = nullptr;
		_j_ptr.clear();
		_w_ptr = nullptr;
	};
	void Species::_copy(const Species& other) {
		_name = other._name;
		_nn_count = other._nn_count;
		_triplet_count = other._triplet_count;
		_count = other._count;
		_t_opt_ptr = other._t_opt_ptr;
		_h_ptr = other._h_ptr;
		_j_ptr = other._j_ptr;
		_w_ptr = other._w_ptr;
	};

	/********************
	Set pointer to the opt time variable
	********************/

	void Species::set_opt_time_ptr(int *t_opt_ptr)
	{
		_t_opt_ptr = t_opt_ptr;
	};

	/********************
	Set h, j ptr
	********************/

	void Species::set_h_ptr(IxnParamTraj *h_ptr) {
		_h_ptr = h_ptr;
	};
	void Species::add_j_ptr(Species* sp, IxnParamTraj *j_ptr) {
		_j_ptr[sp] = j_ptr;
	};
	void Species::set_w_ptr(IxnParamTraj *w_ptr) {
		_w_ptr = w_ptr;
	};

	/********************
	Initialize counts for given other species
	********************/

	void Species::count_nn_for_species(std::vector<Species*> sp_vec) {
		for (auto sp: sp_vec) {
			_nn_count[sp] = 0;
		};
	};
	void Species::count_triplets_for_species(std::vector<std::pair<Species*,Species*>> sp_vec) {
		for (auto sp_pair: sp_vec) {
			_triplet_count[sp_pair.first][sp_pair.second] = 0;
			_triplet_count[sp_pair.second][sp_pair.first] = 0;
		};
	};

	/********************
	Validate setup
	********************/

	void Species::validate_setup() const {
		std::cout << "--- Validate species: " << _name << " ---" << std::endl;
		std::cout << "   NNs: " << std::flush;
		for (auto p: _nn_count) {
			std::cout << p.first->name() << " " << std::flush;
		};
		std::cout << std::endl;
		if (!_t_opt_ptr) {
			std::cerr << "ERROR: no time ptr set" << std::endl;
			exit(EXIT_FAILURE);
		} else {
			std::cout << "   Time ptr is set" << std::endl;
		};
		if (!_h_ptr) {
			std::cerr << "ERROR: no h ptr set" << std::endl;
			exit(EXIT_FAILURE);
		} else {
			std::cout << "   h ptr is set" << std::endl;
		};
	};

	/********************
	Setters/getters
	********************/

	double Species::h() const {
		return _h_ptr->get_at_time(*_t_opt_ptr);
	};
	double Species::j(Species *other) const {
		if (_j_ptr.at(other)) {
			return _j_ptr.at(other)->get_at_time(*_t_opt_ptr);
		} else {
			return 0.;
		};
	};
	double Species::w() const {
		return _w_ptr->get_at_time(*_t_opt_ptr);
	};
	int Species::count() const {
		return _count;
	};
	int Species::nn_count(Species *other) const {
		return _nn_count.at(other);
	};
	int Species::triplet_count(Species* other1, Species* other2) const {
		return _triplet_count.at(other1).at(other2);
	};

	std::string Species::name() const { return _name; };

	/********************
	Increment counts
	********************/

	void Species::count_plus() { _count++; };
	void Species::count_minus() { _count--; };
	void Species::nn_count_plus(Species* other) { _nn_count[other]++; };
	void Species::nn_count_minus(Species* other) { _nn_count[other]--; };
	void Species::triplet_count_plus(Species* other1, Species* other2) { 
		_triplet_count[other1][other2]++;
		_triplet_count[other2][other1]++;
	};
	void Species::triplet_count_minus(Species* other1, Species* other2) { 
		_triplet_count[other1][other2]--;
		_triplet_count[other2][other1]--;
	};

	/********************
	Reset counts
	********************/

	void Species::reset_counts() {
		_count = 0;
		for (auto it=_nn_count.begin(); it!=_nn_count.end(); it++) {
			it->second = 0;
		};
		for (auto it1=_triplet_count.begin(); it1!=_triplet_count.end(); it1++) {
			for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
				it2->second = 0;
			};
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

	/********************
	Constructor
	********************/

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
	Getters
	********************/

	std::string HiddenSpecies::name() const {
		return _name;
	};

};