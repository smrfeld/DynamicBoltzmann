// string
#ifndef STRING_h
#define STRING_h
#include <string>
#endif

// utility for pair
#ifndef PAIR_h
#define PAIR_h
#include <utility>
#endif

// map
#ifndef MAP_h
#define MAP_h
#include <map>
#endif

// Counter
#ifndef COUNTER_h
#define COUNTER_h
#include "counter.hpp"
#endif

/************************************
* Namespace for Gillespie3D
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Doublets, Triplets of Species ptrs
	****************************************/

	class HiddenSpecies;

	struct Species2 {
		Species *sp1;
		Species *sp2;
		Species2(Species* sp1, Species* sp2);
	};
	// Comparator
	bool operator <(const Species2& a, const Species2& b);

	struct Species3 {
		Species *sp1;
		Species *sp2;
		Species *sp3;
		Species3(Species* sp1, Species* sp2, Species *sp3);
	};
	// Comparator
	bool operator <(const Species3& a, const Species3& b);

	struct Species4 {
		Species *sp1;
		Species *sp2;
		Species *sp3;
		Species *sp4;
		Species4(Species* sp1, Species* sp2, Species *sp3, Species *sp4);
	};
	// Comparator
	bool operator <(const Species4& a, const Species4& b);


	struct SpeciesVH {
		Species *sp_visible;
		HiddenSpecies *sp_hidden;

		SpeciesVH(Species* sp_visible, HiddenSpecies *sp_hidden);
	};
	// Comparator
	bool operator <(const SpeciesVH& a, const SpeciesVH& b);

	/****************************************
	Species
	****************************************/
	
	class IxnParam;

	class Species {

	private:

		// Name
		std::string _name;

		// Counters involved
		Counter* _count;
		std::map<Species*, Counter*> _nn_count;
		std::map<Species2, Counter*> _triplet_count;
		std::map<Species3, Counter*> _quartic_count;

		// Pointers to the interaction params
		std::vector<IxnParam*> _h_ptrs;
		std::map<Species*,std::vector<IxnParam*> > _j_ptrs;
		std::map<Species2,std::vector<IxnParam*>> _k_ptrs;
		std::map<Species3,std::vector<IxnParam*>> _q_ptrs;

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const Species& other);

	public:

		/********************
		Constructor
		********************/

		Species(std::string name);
		Species(const Species& other);
		Species(Species&& other);
		Species& operator=(const Species& other);
		Species& operator=(Species&& other);
		~Species();

		/********************
		Set counters
		********************/

		void set_counter(Counter *ctr);
		void add_nn_counter(Species *sp, Counter *ctr);
		void add_triplet_counter(Species *sp1, Species *sp2, Counter *ctr);
		void add_quartic_counter(Species *sp1, Species *sp2, Species *sp3, Counter *ctr);

		/********************
		Set ptr for h,J,K
		********************/

		void add_h_ptr(IxnParam *h_ptr);
		void add_j_ptr(Species* sp, IxnParam *j_ptr);
		void add_k_ptr(Species* sp1, Species* sp2, IxnParam *k_ptr);
		void add_q_ptr(Species* sp1, Species* sp2, Species *sp3, IxnParam *q_ptr);

		/********************
		Ixn params
		********************/

		double h() const;
		double j(Species* other) const;
		double k(Species* other1, Species *other2) const;
		double q(Species* other1, Species *other2, Species *other3) const;

		/********************
		Counts
		********************/

		double count() const;
		double nn_count(Species* other) const;
		double triplet_count(Species* other1, Species *other2) const;
		double quartic_count(Species* other1, Species *other2, Species *other3) const;

		/********************
		Name
		********************/

		std::string name() const;

		/********************
		Increment counts
		********************/

		void count_increment(double inc);
		void nn_count_increment(Species* other, double inc);
		void triplet_count_increment(Species* other1, Species* other2, double inc);
		void quartic_count_increment(Species* other1, Species* other2, Species *other3, double inc);

		/********************
		Reset counts
		********************/

		void reset_counts();
	};
	// Comparator
	bool operator <(const Species& a, const Species& b);

	/****************************************
	HiddenSpecies
	****************************************/
	
	class HiddenSpecies {

	private:

		// Name
		std::string _name;

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const HiddenSpecies& other);

	public:

		/********************
		Constructor
		********************/

		HiddenSpecies(std::string name);
		HiddenSpecies(const HiddenSpecies& other);
		HiddenSpecies(HiddenSpecies&& other);
		HiddenSpecies& operator=(const HiddenSpecies& other);
		HiddenSpecies& operator=(HiddenSpecies&& other);
		~HiddenSpecies();

		/********************
		Name
		********************/

		std::string name() const;
	};
	// Comparator
	bool operator <(const HiddenSpecies& a, const HiddenSpecies& b);

};

