#include "general.hpp"
#include <iomanip>
#include <sstream>

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/********************
	Zero pad a string
	********************/

	std::string pad_str(int i, int n_zeros) {
		std::stringstream fname;
		fname << std::setfill('0') << std::setw(n_zeros) << i;
		return fname.str();
	};

	/********************
	Random numbers
	********************/

	double randD(double dMin, double dMax)
	{
	    return dMin + ((double)rand() / RAND_MAX) * (dMax - dMin);
	};

	int randI(int iMin, int iMax)
	{
		return iMin + rand() % (iMax - iMin + 1);
	};

};