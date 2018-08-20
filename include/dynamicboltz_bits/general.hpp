#ifndef STRING_h
#define STRING_h
#include <string>
#endif

#ifndef VECTOR_h
#define VECTOR_h
#include <vector>
#endif

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	General
	****************************************/

	/********************
	Zero pad a string
	********************/

	std::string pad_str(int i, int n_zeros);
	
	/********************
	Random numbers
	********************/

	double randD(double dMin, double dMax); // inclusive
	int randI(int iMin, int iMax); // inclusive

	/********************
	Pointer deletions
	********************/

	template <typename T> void safeDel(T*& p) { 
		if (p != nullptr) {
			delete p; p=nullptr;
		};
	};
	template <typename T> void safeDelArr(T*& p) { 
		if (p != nullptr) {
			delete[] p; p=nullptr;
		};
	};

	/********************
	Sample a vector of propensities (cumulative probabilities)		
	********************/

	int sample_prop_vec(std::vector<double> &props);

};