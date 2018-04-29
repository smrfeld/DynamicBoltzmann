#ifndef LATTICE_h
#define LATTICE_h
#include "lattice.hpp"
#endif

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/************************************
	Hidden unit
	************************************/

	class HiddenUnit
	{
	private:

		// Species possible
		std::vector<HiddenSpecies*> _sp_possible;

		// Connections
		std::vector<ConnectionVH*> _conn;

		// Biases
		std::map<HiddenSpecies*, std::vector<IxnParam*> > _bias;

		// Probs
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