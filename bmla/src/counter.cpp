#include <iostream>
#include "species.hpp" // includes counter.hpp
#include <algorithm>
#include <iomanip>
#include <fstream>

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
		_report_during_sampling = false;

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
		_counts_stored.clear();
		_counts_stored_averaged.clear();
		_report_during_sampling = false;
	};
	void Counter::_copy(const Counter& other) {
		_type = other._type;
		_binary = other._binary;
		_sp1 = other._sp1;
		_sp2 = other._sp2;
		_sp3 = other._sp3;
		_sp4 = other._sp4;
		_count = other._count;
		_counts_stored = other._counts_stored;
		_counts_stored_averaged = other._counts_stored_averaged;
		_report_during_sampling = other._report_during_sampling;
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
	bool Counter::is_counting_species(std::string s) const {
		if (_sp1) {
			if (_sp1->name() == s) {
				return true;
			};
		};
		return false;
	};
	bool Counter::is_counting_species(std::string s1, std::string s2) const {
		if (_sp1 && _sp2) {
			std::vector<std::string> sps;
			sps.push_back(_sp1->name());
			sps.push_back(_sp2->name());
			auto it1 = std::find(sps.begin(),sps.end(),s1);
			auto it2 = std::find(sps.begin(),sps.end(),s2);;
			if (it1 != sps.end() && it2 != sps.end()) {
				return true;
			};
		};
		return false;
	};
	bool Counter::is_counting_species(std::string s1, std::string s2, std::string s3) const {
		if (_sp1 && _sp2 && _sp3) {
			std::vector<std::string> sps;
			sps.push_back(_sp1->name());
			sps.push_back(_sp2->name());
			sps.push_back(_sp3->name());
			auto it1 = std::find(sps.begin(),sps.end(),s1);
			auto it2 = std::find(sps.begin(),sps.end(),s2);
			auto it3 = std::find(sps.begin(),sps.end(),s3);
			if (it1 != sps.end() && it2 != sps.end() && it3 != sps.end()) {
				return true;
			};
		};
		return false;
	};
	bool Counter::is_counting_species(std::string s1, std::string s2, std::string s3, std::string s4) const {
		if (_sp1 && _sp2 && _sp3 && _sp4) {
			std::vector<std::string> sps;
			sps.push_back(_sp1->name());
			sps.push_back(_sp2->name());
			sps.push_back(_sp3->name());
			sps.push_back(_sp4->name());
			auto it1 = std::find(sps.begin(),sps.end(),s1);
			auto it2 = std::find(sps.begin(),sps.end(),s2);
			auto it3 = std::find(sps.begin(),sps.end(),s3);
			auto it4 = std::find(sps.begin(),sps.end(),s4);
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
	void Counter::reset_count() {
		_count = 0.0;
	};

	// Whether to report this moment during sampling
	void Counter::set_report_during_sampling(bool flag) {
		_report_during_sampling = flag;
	};
	bool Counter::report_during_sampling() const {
		return _report_during_sampling;
	};


	// Committ current count to storage
	void Counter::storage_clear() {
		_counts_stored.clear();
	};
	void Counter::storage_committ_current_count() {
		_counts_stored.push_back(_count);
	};
	double Counter::storage_get_ave_count() const {
		if (_counts_stored.size() == 0) { return 0.; };
		double t=0.0;
		for (auto x: _counts_stored) {
			t += x;
		};
		return t/_counts_stored.size();
	};

	// Average the storage
	void Counter::storage_averaged_committ_current_traj(int average_size) {
		if (_counts_stored_averaged.size() != _counts_stored.size()) {
			// Replace
			_counts_stored_averaged.clear();
			for (int i=0; i<_counts_stored.size(); i++) {
				_counts_stored_averaged.push_back( _counts_stored[i]/average_size);
			};
		} else {
			for (int i=0; i<_counts_stored.size(); i++) {
				_counts_stored_averaged[i] += _counts_stored[i]/average_size;
			};
		};
	};
	void Counter::storage_averaged_clear() {
		_counts_stored_averaged.clear();
	};
	void Counter::storage_averaged_print() const {
		std::cout << "--- Counter: ";
		if (_sp1) {
			std::cout << _sp1->name() << " ";
		};
		if (_sp2) {
			std::cout << _sp2->name() << " ";
		};
		if (_sp3) {
			std::cout << _sp3->name() << " ";
		};
		if (_sp4) {
			std::cout << _sp4->name() << " ";
		};
		std::cout << "---" << std::endl;
		if (_counts_stored_averaged.size() > 0) {
			std::cout << "ave: " << storage_get_ave_count() << " final: " << _counts_stored_averaged.back() << std::endl;
		};
	};
	void Counter::storage_averaged_write(std::ofstream &f) {
		if (_type == COUNT) {
			f << "count " << _sp1->name();
		} else if (_type == NN) {
			f << "NN " << _sp1->name() << " " << _sp2->name();
		} else if (_type == TRIPLET) {
			f << "triplet " << _sp1->name() << " " << _sp2->name() << " " << _sp3->name();
		} else if (_type == QUARTIC) {
			f << "quartic " << _sp1->name() << " " << _sp2->name() << " " << _sp3->name() << " " << _sp4->name();
		};
		for (auto x: _counts_stored_averaged) {
			f << " " << x;
		};
		f << std::endl;
	};
};