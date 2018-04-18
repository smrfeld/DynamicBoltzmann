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

/************************************
* Namespace for Gillespie3D
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Forward declare
	****************************************/

	class IxnParam;
	class HiddenUnit;

	/****************************************
	Species
	****************************************/
	
	class Species {

	private:

		// Name
		std::string _name;

		// Counts
		std::map<Species*, std::map<Species*, double>> _triplet_count;
		std::map<Species*,double> _nn_count;
		double _count;

		// Pointers to the interaction params
		IxnParam *_h_ptr;
		std::map<Species*,IxnParam*> _j_ptr;
		std::map<Species*,std::map<Species*,IxnParam*>> _k_ptr;
		IxnParam *_w_ptr;

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const Species& other);

	public:

		// Constructor
		Species(std::string name);
		Species(const Species& other);
		Species(Species&& other);
		Species& operator=(const Species& other);
		Species& operator=(Species&& other);
		~Species();

		// Initialize counts by informing this species of the existence of all others
		void init_counts(std::list<Species>& sp_list);

		// Set h, j ptr
		void set_h_ptr(IxnParam *h_ptr);
		void add_j_ptr(Species* sp, IxnParam *j_ptr);
		void set_w_ptr(IxnParam *w_ptr);
		void add_k_ptr(Species* sp1, Species* sp2, IxnParam *k_ptr);

		// Setters/getters
		double h() const;
		double j(Species* other) const;
		double w() const;
		double k(Species* other1, Species *other2) const;
		double count() const;
		double nn_count(Species* other) const;
		double triplet_count(Species* other1, Species *other2) const;
		std::string name() const;

		// Increment counts
		void count_increment(double inc);
		void nn_count_increment(Species* other, double inc);
		void triplet_count_increment(Species* other1, Species* other2, double inc);

		// Reset counts
		void reset_counts();
	};
	// Comparator
	bool operator <(const Species& a, const Species& b);
};

