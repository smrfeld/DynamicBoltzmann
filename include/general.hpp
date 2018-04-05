// string
#ifndef STRING_h
#define STRING_h
#include <string>
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
};