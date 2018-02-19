/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Lattice
	****************************************/

	class Lattice
	{
	private:

		// Range in 1d
		double _nu_min;
		double _nu_max;
		
		// Number increments
		int _n_nu;

		// Increment size
		double _dnu;

		// Grid to store solution
		std::map<double,double> _map;



	};

};