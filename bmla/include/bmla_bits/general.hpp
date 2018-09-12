#ifndef STRING_H
#define STRING_H
#include <string>
#endif

#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

/************************************
* Namespace for bmla
************************************/

namespace bmla {

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

	double randD(double dMin, double dMax);
	int randI(int iMin, int iMax);

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