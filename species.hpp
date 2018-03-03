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

		// Current time in the optimization
		int *_t_opt_ptr;

		// Pointers to the interaction params
		IxnParam *_h_ptr;
		std::map<Species*,IxnParam*> _j_ptr;

	public:

		// Constructor
		Species(std::string name);

		// Set pointer to the opt time variable
		void set_opt_time_ptr(int *t_opt_ptr);

		// Set h, j ptr
		void set_h_ptr(IxnParam *h_ptr);
		void add_j_ptr(Species* sp, IxnParam *j_ptr);

		// Validate setup
		void validate_setup() const;

		// Setters/getters
		double h() const;
		double j(Species* other) const;
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

