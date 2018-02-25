#include <iostream>
#include "basis_func_nd.hpp"

using namespace DynamicBoltzmann;

int main() {

	/********************
	Test Grid
	********************/

	/*

	// Make some dimensions
	int n_dim = 3;
	std::vector<Dim> dims;
	dims.push_back(Dim("x",0,1,10));
	dims.push_back(Dim("y",0,1,5));
	dims.push_back(Dim("z",0,1,7));

	// Make grid
	GridND g(n_dim,dims);

	// Access
	int *idxs = new int[n_dim];
	idxs[0] = 0;
	idxs[1] = 2;
	idxs[2] = 1;
	std::cout << g.get(idxs) << std::endl;

	idxs = g.get_idxs(15);
	for (int i=0; i<n_dim; i++) {std::cout << idxs[i] << std::endl;};

	delete[] idxs;
	
	*/

	/****************************************
	Test BasisFuncND
	****************************************/

	/********************
	Test 2D
	********************/

	// Make some dimensions
	std::vector<Dim> dims2d;
	dims2d.push_back(Dim("x",0,4,10));
	dims2d.push_back(Dim("y",0,8,10));

	// Make basis func
	BasisFuncND bf2d(2,dims2d);

	// Fill with sin
	bf2d.test_fill_2d();

	// Evaluate at point
	double x[2] = {2.2,4.5};
	std::cout << bf2d.get_val(x) << " theory " << bf2d.test_func_2d(x[0],x[1]) << std::endl;

	/********************
	Test 2D
	********************/

	// Make some dimensions
	std::vector<Dim> dims3d;
	dims3d.push_back(Dim("x",0,4,10));
	dims3d.push_back(Dim("y",0,8,10));
	dims3d.push_back(Dim("z",-5,2,33));

	// Make basis func
	BasisFuncND bf3d(3,dims3d);

	// Fill with sin
	bf3d.test_fill_3d();

	// Evaluate at point
	double y[3] = {2.2,4.5,2};
	std::cout << bf3d.get_val(y) << " theory " << bf3d.test_func_3d(y[0],y[1],y[2]) << std::endl;

	return 0;
};