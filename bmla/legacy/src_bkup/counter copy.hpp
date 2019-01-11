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
	Counter
	****************************************/
	
	class Species; 

	enum CounterType { COUNT, NN, TRIPLET, QUARTIC };

	class Counter {

	private:

		// Type
		CounterType _type;

		// Binary or probabilistic
		bool _binary;

		// Species involved
		Species *_sp1,*_sp2,*_sp3,*_sp4;

		// Counts as double
		// If binary, can be cast!
		double _count;

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const Counter& other);

	public:

		// Constructor
		Counter(Species *sp1, bool binary=true);
		Counter(Species *sp1, Species *sp2, bool binary=true);
		Counter(Species *sp1, Species *sp2, Species *sp3, bool binary=true);
		Counter(Species *sp1, Species *sp2, Species *sp3, Species *sp4, bool binary=true);
		Counter(const Counter& other);
		Counter(Counter&& other);
		Counter& operator=(const Counter& other);
		Counter& operator=(Counter&& other);
		~Counter();

		// Switch binary/prob
		void set_binary(bool flag);

		// Check binary
		bool is_binary() const;

		// Check type
		bool is_type(CounterType counter_type) const;

		// Check species
		bool counts_species(std::string s) const;
		bool counts_species(std::string s1, std::string s2) const;
		bool counts_species(std::string s1, std::string s2, std::string s3) const;
		bool counts_species(std::string s1, std::string s2, std::string s3, std::string s4) const;

		// Get count
		double get_count() const;

		// Increment count
		void increment(double inc);

		// Reset counts
		void reset_counts();
	};

	/****************************************
	Container structure to find counters easily, etc.
	****************************************/

	class CounterContainer {
	private:

		// List of counters
		std::list<Counter> _counters;

		// Maps to find by species
		std::map<Species*, Counter*> _map1;
		std::map<Species*, std::map<Species*, Counter*>> _map2;
		std::map<Species*, std::map<Species*, std::map<Species*, Counter*>>> _map3;
		std::map<Species*, std::map<Species*, std::map<Species*, std::map<Species*, Counter*>>>> _map4;

		// Add to maps
		void _add_to_maps(std::set<Species*> sp_set, Counter *counter);

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const CounterContainer& other);

	public:

		// Constructor
		CounterContainer();
		CounterContainer(const CounterContainer& other);
		CounterContainer(CounterContainer&& other);
		CounterContainer& operator=(const CounterContainer& other);
		CounterContainer& operator=(CounterContainer&& other);
		~CounterContainer();

		/********************
		Get a counter by species ptr(s)
		********************/

		// Note: these really shouldnt be used, but needed by IxnParams!
		Counter* get_counter(Species *sp) const;
		Counter* get_counter(Species *sp1, Species *sp2) const;
		Counter* get_counter(Species *sp1, Species *sp2, Species *sp3) const;
		Counter* get_counter(Species *sp1, Species *sp2, Species *sp3, Species *sp4) const;

		/********************
		Get counts by species ptr(s)
		********************/

		double get_count(Species *sp) const;
		double get_count(Species *sp1, Species *sp2) const;
		double get_count(Species *sp1, Species *sp2, Species *sp3) const;
		double get_count(Species *sp1, Species *sp2, Species *sp3, Species *sp4) const;
		int get_count_binary(Species *sp) const;
		int get_count_binary(Species *sp1, Species *sp2) const;
		int get_count_binary(Species *sp1, Species *sp2, Species *sp3) const;
		int get_count_binary(Species *sp1, Species *sp2, Species *sp3, Species *sp4) const;

		/********************
		Increment count(s) by species
		********************/

		void count_increment(Species *sp, double inc);
		void count_increment(Species *sp1, Species *sp2, double inc);
		void count_increment(Species *sp1, Species *sp2, Species *sp3, double inc);
		void count_increment(Species *sp1, Species *sp2, Species *sp3, Species *sp4, double inc);
		void count_plus(Species *sp);
		void count_plus(Species *sp1, Species *sp2);
		void count_plus(Species *sp1, Species *sp2, Species *sp3);
		void count_plus(Species *sp1, Species *sp2, Species *sp3, Species *sp4);
		void count_minus(Species *sp);
		void count_minus(Species *sp1, Species *sp2);
		void count_minus(Species *sp1, Species *sp2, Species *sp3);
		void count_minus(Species *sp1, Species *sp2, Species *sp3, Species *sp4);

		/********************
		Make a counter
		********************/

		void add_counter(Species *sp, bool binary=true);
		void add_counter(Species *sp1, Species *sp2, bool binary=true);
		void add_counter(Species *sp1, Species *sp2, Species *sp3, bool binary=true);
		void add_counter(Species *sp1, Species *sp2, Species *sp3, Species *sp4, bool binary=true);
	};
};

