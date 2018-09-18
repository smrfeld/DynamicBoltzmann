#ifndef STRING_h
#define STRING_h
#include <string>
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	General
	****************************************/

	/********************
	Sign function
	********************/

	int sgn(double val);

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
};