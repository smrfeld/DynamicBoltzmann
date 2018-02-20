#include "timeseries.hpp"
#include "lattice.hpp"

using namespace DynamicBoltzmann;

int main() {

	// Timeseries t(0.0,1.0,10);
	// std::cout << t << std::endl;

	Lattice l(0.0,1.0,11);
	l.test_fill();
	std::cout << l << std::endl;

	l.write_to_file("l1.txt");
	l.write_to_file("l2.txt",1000);
	l.write_deriv_to_file("l3.txt");
	l.write_deriv_to_file("l4.txt",1000);

	return 0;
};