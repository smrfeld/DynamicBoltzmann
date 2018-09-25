# Dynamic Boltzmann Machine for Lattice Chemical Kinetics

## Dependencies

d-Cubic library for doing cubic interpolation and derivatives of the interpolation in d dimensions.

Get it from [here](https://github.com/smrfeld/d-Cubic).

## Installation

There now is a convenient makefile:
```
make
make install
```
The default install locations are `/usr/local/include` and `/usr/local/lib`.

If not installing, make sure to move the library to somewhere it will be found, or at worst
```
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/absolute/path/to/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/Users/oernst/Research/cnl/dynamic_boltzmann_cpp/lib
```
In this case, you will also need to specify using `-L` the location of the library at linking time (`-L` is needed for finding the library at linking time; `DYLD_LIBRARY_PATH` is at runtime).

## Usage

To access the headers, simply use: `include <dynamicboltz>`. 

Link your program using:
```
g++ -std=c++14 -O3 -ldynamicboltz -ldcubic myprogram.cpp -o myprogram.o
```
(this assumes you have install `-ldcubic`).

## Namespace

The namespace is `dblz` ("DynamicBoltzmann").

# BMLA (LEGACY)

NOTE: THIS NEEDS AN UPDATE!

To build the library called `bmla` - navigate to the dir `bmla`, then as above! A convenient headir is `/usr/local/include/bmla`.

## Compilation

Link your program using:
```
g++ -std=c++14 -O3 -lbmla myprogram.cpp -o myprogram.o
```

## Namespace

The namespace is `bmla`.