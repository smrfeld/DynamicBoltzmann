# BMLA for initial conditions

## Dependencies

None.

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

To access the headers, simply use: `include <bmla>`. 

Link your program using:
```
g++ -std=c++14 -O3 -lbmla bmla.cpp -o bmla.o
```

## Namespace

The namespace is `bmla`.