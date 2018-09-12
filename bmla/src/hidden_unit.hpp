#ifndef VECTOR_H
#define VECTOR_H
#include <vector> 
#endif

#ifndef MAP_H
#define MAP_H
#include <map> 
#endif

/************************************
* Namespace for bmla
************************************/

namespace bmla {

	// Forward
	class ConnectionVH;
	class HiddenSpecies;
	class IxnParam;

	/************************************
	Hidden unit
	************************************/

	class HiddenUnit
	{
	private:

		// Connections
		std::vector<ConnectionVH*> _conn;

		// Biases
		std::map<HiddenSpecies*, std::vector<IxnParam*> > _bias;

		// Probs
		std::vector<HiddenSpecies*> _sp_possible;
		std::map<HiddenSpecies*, double> _probs;
		double _prob_empty;

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const HiddenUnit& other);

	public:

		/********************
		Constructor
		********************/

		HiddenUnit();
		HiddenUnit(const HiddenUnit& other);
		HiddenUnit(HiddenUnit&& other);
		HiddenUnit& operator=(const HiddenUnit& other);
		HiddenUnit& operator=(HiddenUnit&& other);
		~HiddenUnit();	

		/********************
		Add possible species
		********************/

		void add_hidden_species_possibility(HiddenSpecies* sp);

		/********************
		Add a connection
		********************/

		void add_visible_hidden_conn(ConnectionVH* conn);

		/********************
		Print connections
		********************/

		void print_conns(bool newline) const;

		/********************
		Getters
		********************/

		// Nullptr for prob of empty
		double get_prob(HiddenSpecies *hsp) const;

		/********************
		Add a bias
		********************/

		void add_bias(IxnParam *ip);

		/********************
		Activate
		********************/

		void activate(bool binary); 

		/********************
		Binarize
		********************/

		void binarize();

	};
};