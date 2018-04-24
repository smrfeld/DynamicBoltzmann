# Dynamic Boltzmann Machine for Lattice Chemical Kinetics

## Build

There now is a convenient makefile:
```
make
make install
```
Headers are installed in `/usr/local/include/dynamicboltz` and libraries in `/usr/local/lib`. To access the headers, make sure to `include "dynamicboltz/....hpp"`. Link using:
```
g++ -std=c++14 -O3 -ldynamicboltz myprogram.cpp -o myprogram.o
```
If not installing, make sure to move the library to somewhere it will be found, or at worst
```
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/absolute/path/to/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/Users/oernst/Research/cnl/dynamic_boltzmann_cpp/lib
```
(Note: `-L` is needed for finding the library at linking time; `DYLD_LIBRARY_PATH` is at runtime).

# BMLA

To build the library called `bmla` - navigate to the dir `bmla`, then as above! Headers are in `/usr/local/include/bmla`.