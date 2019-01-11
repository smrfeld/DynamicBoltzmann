#ifndef STRING_H
#define STRING_H
#include <string>
#endif

#ifndef VECTOR_H
#define VECTOR_H
#include <vector> 
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {

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

		// Storage of counts, as desired
		std::vector<double> _counts_stored, _counts_stored_averaged;

		// Whether this counter should be reported during sampling
		bool _report_during_sampling;

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const Counter& other);

	public:

		/********************
		Constructor
		********************/

		Counter(Species *sp1, bool binary=true);
		Counter(Species *sp1, Species *sp2, bool binary=true);
		Counter(Species *sp1, Species *sp2, Species *sp3, bool binary=true);
		Counter(Species *sp1, Species *sp2, Species *sp3, Species *sp4, bool binary=true);
		Counter(const Counter& other);
		Counter(Counter&& other);
		Counter& operator=(const Counter& other);
		Counter& operator=(Counter&& other);
		~Counter();

		/********************
		Switch binary/prob
		********************/

		void set_binary(bool flag);

		/********************
		Check binary
		********************/

		bool is_binary() const;

		/********************
		Check type
		********************/

		bool is_type(CounterType counter_type) const;

		/********************
		Check species
		********************/

		bool is_counting_species(std::string s) const;
		bool is_counting_species(std::string s1, std::string s2) const;
		bool is_counting_species(std::string s1, std::string s2, std::string s3) const;
		bool is_counting_species(std::string s1, std::string s2, std::string s3, std::string s4) const;

		/********************
		Get count
		********************/

		double get_count() const;

		/********************
		Increment count
		********************/

		void increment(double inc);

		/********************
		Reset current count
		********************/

		void reset_count();

		/********************
		Whether to report this moment during sampling
		********************/

		void set_report_during_sampling(bool flag);
		bool report_during_sampling() const;

		/********************
		Committ current count to storage
		********************/

		void storage_clear();
		void storage_committ_current_count();

		/********************
		Average the storage
		********************/

		void storage_averaged_committ_current_traj(int average_size);
		double storage_averaged_get_ave_count() const;
		void storage_averaged_clear();
		void storage_averaged_print() const;
		void storage_averaged_write(std::ofstream &f);
	};

};

