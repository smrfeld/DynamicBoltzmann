# Dynamic Boltzmann Machine for Lattice Chemical Kinetics

## Dependencies

None?

## Installation

There now is a convenient makefile:
```
make
make install
```
Headers are by default installed in `/usr/local/include/dynamicboltz_bits` with a convenient header `/usr/local/include/dynamicboltz` and libraries in `/usr/local/lib`. To access the headers, simply use: `include <dynamicboltz>`. 

If not installing, make sure to move the library to somewhere it will be found, or at worst
```
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/absolute/path/to/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/Users/oernst/Research/cnl/dynamic_boltzmann_cpp/lib
```
(Note: `-L` is needed for finding the library at linking time; `DYLD_LIBRARY_PATH` is at runtime).

## Compilation

Link your program using:
```
g++ -std=c++14 -O3 -ldynamicboltz myprogram.cpp -o myprogram.o
```

## Namespace

The namespace is `dboltz` ("DynamicBoltzmann").

# BMLA

To build the library called `bmla` - navigate to the dir `bmla`, then as above! Headers are in `/usr/local/include/bmla`.