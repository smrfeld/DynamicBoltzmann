#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

#ifndef MAP_H
#define MAP_H
#include <map>
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/************************************
	Hidden unit
	************************************/

	// Forward
	class ConnectionVH;
	class HiddenSpecies;
	class IxnParamTraj;

	class UnitHidden
	{
	private:

		// Dimension and location
		int _layer;
		int _dim;
		int _x,_y,_z;

		// Connections
		std::vector<ConnectionVH*> _conn;

		// Biases
		std::map<HiddenSpecies*, std::vector<IxnParamTraj*> > _bias;

		// Probs
		std::vector<HiddenSpecies*> _sp_possible;
		std::map<HiddenSpecies*, double> _probs;
		double _prob_empty;

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const UnitHidden& other);

	public:

		/********************
		Constructor
		********************/

		UnitHidden();
		UnitHidden(const UnitHidden& other);
		UnitHidden(UnitHidden&& other);
		UnitHidden& operator=(const UnitHidden& other);
		UnitHidden& operator=(UnitHidden&& other);
		~UnitHidden();	

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

		void add_bias(IxnParamTraj *ip);

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