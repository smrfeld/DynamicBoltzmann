#include <iostream>
#include "counter.hpp"

/************************************
* Namespace for Gillespie3D
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Counter
	****************************************/

	Counter::Counter(Species *sp1, bool binary) : Counter(sp1, nullptr, nullptr, nullptr, binary) {
		_type = COUNT;
	};
	Counter::Counter(Species *sp1, Species *sp2, bool binary) : Counter(sp1, sp2, nullptr, nullptr, binary) {
		_type = NN;
	};
	Counter::Counter(Species *sp1, Species *sp2, Species *sp3, bool binary) : Counter(sp1, sp2, sp3, nullptr, binary) {
		_type = TRIPLET;
	};
	Counter::Counter(Species *sp1, Species *sp2, Species *sp3, Species *sp4, bool binary) {
		_sp1 = sp1;
		_sp2 = sp2;
		_sp3 = sp3;
		_sp4 = sp4;
		_binary = binary;
		_type = QUARTIC;

		// Counts
		_count = 0.0;
	};
	Counter::Counter(const Counter& other) {
		_copy(other);
	};
	Counter::Counter(Counter&& other) {
		_copy(other);
		other._reset();
	};
	Counter& Counter::operator=(const Counter& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	Counter& Counter::operator=(Counter&& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
	};
	Counter::~Counter() {
		_clean_up();
	};

	void Counter::_clean_up() {
		// Nothing...
	};
	void Counter::_reset() {
		_binary = true;
		_sp1 = nullptr;
		_sp2 = nullptr;
		_sp3 = nullptr;
		_sp4 = nullptr;
		_count = 0.0;
	};
	void Counter::_copy(const Counter& other) {
		_type = other._type;
		_binary = other._binary;
		_sp1 = other._sp1;
		_sp2 = other._sp2;
		_sp3 = other._sp3;
		_sp4 = other._sp4;
		_count = other._count;
	};

	// Switch binary/prob
	void Counter::set_binary(bool flag) {
		_binary = flag;
	};

	// Check binary
	bool Counter::is_binary() const {
		return _binary;
	};

	// Check type
	bool Counter::is_type(CounterType counter_type) const {
		if (_type == counter_type) {
			return true;
		} else {
			return false;
		};	
	};

	// Check species
	bool Counter::counts_species(std::string s) const {
		if (_sp1)
			if (_sp1->name() == s) {
				return true;
			};
		};
		return false;
	};
	bool Counter::counts_species(std::string s1, std::string s2) const {
		if (_sp1 && _sp2) {
			std::vector<std::string> sps;
			sps.push_back(_sp1->name());
			sps.push_back(_sp2->name());
			auto it1 = sps.find(s1);
			auto it2 = sps.find(s2);
			if (it1 != sps.end() && it2 != sps.end()) {
				return true;
			};
		};
		return false;
	};
	bool Counter::counts_species(std::string s1, std::string s2, std::string s3) const {
		if (_sp1 && _sp2 && _sp3) {
			std::vector<std::string> sps;
			sps.push_back(_sp1->name());
			sps.push_back(_sp2->name());
			sps.push_back(_sp3->name());
			auto it1 = sps.find(s1);
			auto it2 = sps.find(s2);
			auto it3 = sps.find(s3);
			if (it1 != sps.end() && it2 != sps.end() && it3 != sps.end()) {
				return true;
			};
		};
		return false;
	};
	bool Counter::counts_species(std::string s1, std::string s2, std::string s3, std::string s4) const {
		if (_sp1 && _sp2 && _sp3 && _sp4) {
			std::vector<std::string> sps;
			sps.push_back(_sp1->name());
			sps.push_back(_sp2->name());
			sps.push_back(_sp3->name());
			sps.push_back(_sp4->name());
			auto it1 = sps.find(s1);
			auto it2 = sps.find(s2);
			auto it3 = sps.find(s3);
			auto it4 = sps.find(s4);
			if (it1 != sps.end() && it2 != sps.end() && it3 != sps.end() && it4 != sps.end()) {
				return true;
			};
		};
		return false;
	};

	// Get count
	double Counter::get_count() const {
		return _count;
	};

	// Increment count
	void Counter::increment(double inc) {
		_count += inc;
	};

	// Reset counts
	void Counter::reset_counts() {
		_count = 0.0;
	};

	/****************************************
	Container structure to find counters easily, etc.
	****************************************/

	CounterContainer::CounterContainer() {};
	CounterContainer::CounterContainer(const CounterContainer& other) {
		_copy(other);
	};
	CounterContainer::CounterContainer(CounterContainer&& other) {
		_copy(other);
		other._reset();
	};
	CounterContainer& CounterContainer::operator=(const CounterContainer& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	CounterContainer& CounterContainer::operator=(CounterContainer&& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
	};
	CounterContainer::~CounterContainer() {
		_clean_up();
	};

	void CounterContainer::_clean_up() {
		// Nothing...
	};
	void CounterContainer::_reset() {
		_counter.clear();
		_map1.clear();
		_map2.clear();
		_map3.clear();
		_map4.clear();
	};
	void CounterContainer::_copy(const CounterContainer& other) {
		_counter = other._counter;
		_map1 = other._map1;
		_map2 = other._map2;
		_map3 = other._map3;
		_map4 = other._map4;
	};

	/********************
	Get a counter by species ptr(s)
	********************/

	Counter* CounterContainer::get_counter(Species *sp) const {
		auto it = _map1.find(sp);
		if (it != _map1.end()) {
			return &(it->second);
		};
		return nullptr;
	};
	Counter* CounterContainer::get_counter(Species *sp1, Species *sp2) const {
		auto it1 = _map2.find(sp1);
		if (it1 != _map2.end()) {
			auto it2 = it1->second.find(sp2);
			if (it2 != it1->second.end()) {
				return &(it2->second);
			};
		};
		return nullptr;
	};
	Counter* CounterContainer::get_counter(Species *sp1, Species *sp2, Species *sp3) const {
		auto it1 = _map2.find(sp1);
		if (it1 != _map2.end()) {
			auto it2 = it1->second.find(sp2);
			if (it2 != it1->second.end()) {
				auto it3 = it2->second.find(sp3);
				if (it3 != it2->second.end()) {
					return &(it3->second);
				};
			};
		};
		return nullptr;
	};
	Counter* CounterContainer::get_counter(Species *sp1, Species *sp2, Species *sp3, Species *sp4) const {
		auto it1 = _map2.find(sp1);
		if (it1 != _map2.end()) {
			auto it2 = it1->second.find(sp2);
			if (it2 != it1->second.end()) {
				auto it3 = it2->second.find(sp3);
				if (it3 != it2->second.end()) {
					auto it4 = it3->second.find(sp4);
					if (it4 != it3->second.end()) {
						return &(it4->second);
					};
				};
			};
		};
		return nullptr;
	};

	/********************
	Get counts
	********************/

	double CounterContainer::get_count(Species *sp) const {
		auto ptr = get_counter(sp);
		if (ptr) {
			return ptr->get_count();
		} else {
			return 0.0;
		};
	};
	double CounterContainer::get_count(Species *sp1, Species *sp2) const {
		auto ptr = get_counter(sp1,sp2);
		if (ptr) {
			return ptr->get_count();
		} else {
			return 0.0;
		};
	};
	double CounterContainer::get_count(Species *sp1, Species *sp2, Species *sp3) const {
		auto ptr = get_counter(sp1,sp2,sp3);
		if (ptr) {
			return ptr->get_count();
		} else {
			return 0.0;
		};
	};
	double CounterContainer::get_count(Species *sp1, Species *sp2, Species *sp3, Species *sp4) const {
		auto ptr = get_counter(sp1,sp2,sp3,sp4);
		if (ptr) {
			return ptr->get_count();
		} else {
			return 0.0;
		};
	};
	int CounterContainer::get_count_binary(Species *sp) const {
		auto ptr = get_counter(sp);
		if (ptr) {
			return ptr->get_count_binary();
		} else {
			return 0;
		};
	};
	int CounterContainer::get_count_binary(Species *sp1, Species *sp2) const {
		auto ptr = get_counter(sp1,sp2);
		if (ptr) {
			return ptr->get_count_binary();
		} else {
			return 0;
		};
	};
	int CounterContainer::get_count_binary(Species *sp1, Species *sp2, Species *sp3) const {
		auto ptr = get_counter(sp1,sp2,sp3);
		if (ptr) {
			return ptr->get_count_binary();
		} else {
			return 0;
		};
	};
	int CounterContainer::get_count_binary(Species *sp1, Species *sp2, Species *sp3, Species *sp4) const {
		auto ptr = get_counter(sp1,sp2,sp3,sp4);
		if (ptr) {
			return ptr->get_count_binary();
		} else {
			return 0;
		};
	};

	/********************
	Increment count(s) by species
	********************/

	void CounterContainer::count_increment(Species *sp, double inc) {
		auto ptr = get_counter(sp);
		if (ptr) {
			ptr->increment(inc);
		};
	};
	void CounterContainer::count_increment(Species *sp1, Species *sp2, double inc) {
		auto ptr = get_counter(sp1,sp2);
		if (ptr) {
			ptr->increment(inc);
		};
	};
	void CounterContainer::count_increment(Species *sp1, Species *sp2, Species *sp3, double inc) {
		auto ptr = get_counter(sp1,sp2,sp3);
		if (ptr) {
			ptr->increment(inc);
		};
	};
	void CounterContainer::count_increment(Species *sp1, Species *sp2, Species *sp3, Species *sp4, double inc) {
		auto ptr = get_counter(sp1,sp2,sp3,sp4);
		if (ptr) {
			ptr->increment(inc);
		};
	};
	void CounterContainer::count_plus(Species *sp) {
		auto ptr = get_counter(sp);
		if (ptr) {
			ptr->plus();
		};
	};
	void CounterContainer::count_plus(Species *sp1, Species *sp2) {
		auto ptr = get_counter(sp1,sp2);
		if (ptr) {
			ptr->plus();
		};
	};
	void CounterContainer::count_plus(Species *sp1, Species *sp2, Species *sp3) {
		auto ptr = get_counter(sp1,sp2,sp3);
		if (ptr) {
			ptr->plus();
		};
	};
	void CounterContainer::count_plus(Species *sp1, Species *sp2, Species *sp3, Species *sp4) {
		auto ptr = get_counter(sp1,sp2,sp3,sp4);
		if (ptr) {
			ptr->plus();
		};
	};
	void CounterContainer::count_minus(Species *sp) {
		auto ptr = get_counter(sp);
		if (ptr) {
			ptr->minus();
		};
	};
	void CounterContainer::count_minus(Species *sp1, Species *sp2) {
		auto ptr = get_counter(sp1,sp2);
		if (ptr) {
			ptr->minus();
		};
	};
	void CounterContainer::count_minus(Species *sp1, Species *sp2, Species *sp3) {
		auto ptr = get_counter(sp1,sp2,sp3);
		if (ptr) {
			ptr->minus();
		};
	};
	void CounterContainer::count_minus(Species *sp1, Species *sp2, Species *sp3, Species *sp4) {
		auto ptr = get_counter(sp1,sp2,sp3,sp4);
		if (ptr) {
			ptr->minus();
		};
	};

	/********************
	Make a counter
	********************/

	Counter* CounterContainer::add_counter(Species *sp, bool binary) {
		_counters.push_back(Counter(sp,binary));
		std::set<Species*> sp_set; // list of unique species
		sp_set.push_back(sp);
		_add_to_maps(sp_set, &_counters.back());
		return &_counters.back();
	};
	Counter* CounterContainer::add_counter(Species *sp1, Species *sp2, bool binary) {
		_counters.push_back(Counter(sp1,sp2,binary));
		std::set<Species*> sp_set; // list of unique species
		sp_set.push_back(sp1);
		sp_set.push_back(sp2);
		_add_to_maps(sp_set, &_counters.back());
		return &_counters.back();
	};
	Counter* CounterContainer::add_counter(Species *sp1, Species *sp2, Species *sp3, bool binary) {	
		_counters.push_back(Counter(sp1,sp2,sp3,binary));
		std::set<Species*> sp_set; // list of unique species
		sp_set.push_back(sp1);
		sp_set.push_back(sp2);		
		sp_set.push_back(sp3);
		_add_to_maps(sp_set, &_counters.back());
		return &_counters.back();
	};
	Counter* CounterContainer::add_counter(Species *sp1, Species *sp2, Species *sp3, Species *sp4, bool binary) {
		_counters.push_back(Counter(sp1,sp2,sp3,sp4,binary));
		std::set<Species*> sp_set; // list of unique species
		sp_set.push_back(sp1);
		sp_set.push_back(sp2);		
		sp_set.push_back(sp3);
		sp_set.push_back(sp4);
		_add_to_maps(sp_set, &_counters.back());
		return &_counters.back();
	};



	/********************
	Add to maps
	********************/

	void CounterContainer::_add_to_maps(std::set<Species*> sp_set, Counter* counter) {
		auto it1 = sp_set.begin();
		while (it1 != sp_set.end()) {
			if (sp_set.size() == 1)
			{
				// add
				_map1[*it1] = counter;
			} else {
				auto it2 = it1;
				it2++;
				while (it2 != sp_set.end()) {
					if (sp_set.size() == 2)
					{
						// add
						_map2[*it1][*it2] = counter;
					} else {
						auto it3 = it2;
						it3++;
						while (it3 != sp_set.end()) {
							if (sp_set.size() == 3) 
							{
								// add
								_map3[*it1][*it2][*it3] = counter;
							} else {
								auto it4 = it3;
								it4++;
								while (it4 != sp_set.end()) {
									if (sp_set.size() == 4)
									{
										// add
										_map4[*it1][*it2][*it3][*it4] = counter;
									} else {
										std::cerr << "Error! Container structure failed" << std::endl;
										exit(EXIT_FAILURE);
									};
									// inc
									it4++;
								};
							};
							// inc
							it3++;
						};
					};
					// inc
					it2++;
				};
			};
			// inc
			it1++;
		};
	};


};