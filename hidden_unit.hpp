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

		// Number of sites I am connected to
		int _n_conn;

		// Sites I am connected to
		std::vector<Site*> _conn;

		// Species that I love
		Species *_sp;

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

		HiddenUnit(std::vector<Site*> conn, Species *sp);
		HiddenUnit(const HiddenUnit& other);
		HiddenUnit(HiddenUnit&& other);
		HiddenUnit& operator=(const HiddenUnit& other);
		HiddenUnit& operator=(HiddenUnit&& other);
		~HiddenUnit();	

		/********************
		Getters
		********************/

		double get() const;

		/********************
		Activate
		********************/

		void activate(bool binary); 

	};
};