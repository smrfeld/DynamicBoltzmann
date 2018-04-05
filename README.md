# Dynamic Boltzmann Machine for Lattice Chemical Kinetics

## Build

There now is a convenient makefile:
```
make
```
The library will be in `lib` and the necessary headers in `include`. Link against it using
```
g++ -std=c++14 -O3 -L./lib -ldynamicboltz myprogram.cpp -o myprogram.o
```
Make sure to move the library to somewhere it will be found, or at worst
```
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/absolute/path/to/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/Users/oernst/Research/cnl/dynamic_boltzmann_cpp/lib
```
(Note: `-L` is needed for finding the library at linking time; `DYLD_LIBRARY_PATH` is at runtime).

## LEGACY build instructions

To build the library:
```
g++ -std=c++14 -O3 -c -fpic basis_func.cpp dynamic_boltzmann.cpp general.cpp grid.cpp ixn_param_traj.cpp lattice.cpp species.cpp var_term_traj.cpp hidden_unit.cpp
g++ -std=c++14 -O3 -shared -o libdynamicboltz.so basis_func.o dynamic_boltzmann.o general.o grid.o ixn_param_traj.o lattice.o species.o var_term_traj.o hidden_unit.o
```

Link against it using
```
g++ -std=c++14 -O3 -L./ -ldynamicboltz myprogram.cpp -o myprogram.o
```

# BMLA

To build the library called `bmla` - navigate to the dir `bmla`, then:
```
g++ -std=c++14 -O3 -c -fpic bmla.cpp ixn_param.cpp species.cpp
g++ -std=c++14 -O3 -shared -o libbmla.so bmla.o ../general.o species.o ixn_param.o ../lattice.o ../hidden_unit.o
```

Link against it using
```
g++ -std=c++14 -O3 -L./ -lbmla myprogram.cpp -o myprogram.o
```