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
		std::map<Species*,int> _nn_count;
		int _count;

		// Pointers to the interaction params
		IxnParam *_h_ptr;
		std::map<Species*,IxnParam*> _j_ptr;
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

		// Set h, j ptr
		void set_h_ptr(IxnParam *h_ptr);
		void add_j_ptr(Species* sp, IxnParam *j_ptr);
		void set_w_ptr(IxnParam *w_ptr);

		// Setters/getters
		double h() const;
		double j(Species* other) const;
		double w() const;
		int count() const;
		int nn_count(Species* other) const;
		std::string name() const;

		// Increment counts
		void count_plus();
		void count_minus();
		void nn_count_plus(Species* other);
		void nn_count_minus(Species* other);

		// Reset counts
		void reset_counts();
	};
	// Comparator
	bool operator <(const Species& a, const Species& b);
};
