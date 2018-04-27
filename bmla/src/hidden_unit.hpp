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

	// Forward declare a site
	class Site;

	class HiddenUnit
	{
	private:

		// Connections
		std::vector<ConnectionVH*> _conn;

		// Biases
		std::vector<IxnParam*> _bias;

		// Value
		double _val;

		// Activation function
		double _sigma(double x) const;

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const HiddenUnit& other);

	public:

		/********************
		Constructor
		********************/

		HiddenUnit(std::vector<IxnParam*> bias);
		HiddenUnit(const HiddenUnit& other);
		HiddenUnit(HiddenUnit&& other);
		HiddenUnit& operator=(const HiddenUnit& other);
		HiddenUnit& operator=(HiddenUnit&& other);
		~HiddenUnit();	

		/********************
		Add a connection
		********************/

		void add_connection(ConnectionVH* conn);

		/********************
		Print connections
		********************/

		void print_conns(bool newline) const;

		/********************
		Getters
		********************/

		double get() const;

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