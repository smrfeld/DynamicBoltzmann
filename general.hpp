
/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	General
	****************************************/

	/********************
	Random numbers
	********************/

	// Random numbers
	inline double randD(double dMin, double dMax)
	{
	    return dMin + ((double)rand() / RAND_MAX) * (dMax - dMin);
	};
	inline int randI(int iMin, int iMax)
	{
		return iMin + rand() % (iMax - iMin + 1);
	};

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