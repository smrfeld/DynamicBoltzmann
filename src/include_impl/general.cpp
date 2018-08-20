#include "../../include/dynamicboltz_bits/general.hpp"

#include <iomanip>
#include <sstream>

/************************************
* Namespace for dboltz
************************************/

namespace dboltz {

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

	/********************
	Sample a vector of propensities (cumulative probabilities)		
	********************/

	int sample_prop_vec(std::vector<double> &props) {
		// Sample RV
		double r = randD(0.0,props.back());

		// Find interval
		for (int i=0; i<props.size()-1; i++) {
			if (props[i] <= r && r <= props[i+1]) {
				return i;
			};
		};
		return 0; // never get here
	};

};