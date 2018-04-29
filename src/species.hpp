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

	class IxnParamTraj;
	class HiddenUnit;

	/****************************************
	Species
	****************************************/
	
	class Species {

	private:

		// Name
		std::string _name;

		// Counts
		std::map<Species*,std::map<Species*,int>> _triplet_count;
		std::map<Species*,int> _nn_count;
		int _count;

		// Current time in the optimization
		int *_t_opt_ptr;

		// Pointers to the interaction params
		std::vector<IxnParamTraj*> _h_ptrs;
		std::map<Species*,std::vector<IxnParamTraj*> > _j_ptrs;

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
		Set optimization time pointer
		********************/

		void set_opt_time_ptr(int *t_opt_ptr);

		/********************
		Set h, J, W ptr
		********************/

		void add_h_ptr(IxnParamTraj *h_ptr);
		void add_j_ptr(Species* sp, IxnParamTraj *j_ptr);

		/********************
		Initialize counts for a given other species
		********************/

		void count_nn_for_species(std::vector<Species*> sp_vec);
		void count_triplets_for_species(std::vector<std::pair<Species*,Species*>> sp_vec);

		/********************
		Validate setup
		********************/

		void validate_setup() const;

		/********************
		Getters
		********************/

		double h() const;
		double j(Species* other) const;
		int count() const;
		int nn_count(Species* other) const;
		int triplet_count(Species* other1, Species* other2) const;
		std::string name() const;

		/********************
		Increment counts
		********************/

		void count_plus();
		void count_minus();
		void nn_count_plus(Species* other);
		void nn_count_minus(Species* other);
		void triplet_count_plus(Species* other1, Species* other2);
		void triplet_count_minus(Species* other1, Species* other2);

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
		Getters
		********************/

		std::string name() const;
	};

};




