#ifndef TIMESERIES_h
#define TIMESERIES_h
#include "timeseries.hpp"
#endif

using namespace DynamicBoltzmann;

int main() {

	Timeseries t(0.0,1.0,10);
	t.p();

	Timeseries s(t);
	s.p();

	Timeseries r = t;
	r.p();
	
	return 0;
};