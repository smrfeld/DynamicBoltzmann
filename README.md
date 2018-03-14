# Dynamic Boltzmann Machine for Lattice Chemical Kinetics

## Build

To build the library:
```
g++ -std=c++14 -O3 -c -fpic basis_func.cpp dynamic_boltzmann.cpp general.cpp grid.cpp ixn_param_traj.cpp lattice.cpp species.cpp var_term_traj.cpp
g++ -std=c++14 -O3 -shared -o libdynamicboltz.so basis_func.o dynamic_boltzmann.o general.o grid.o ixn_param_traj.o lattice.o species.o var_term_traj.o
```

Link against it using
```
g++ -std=c++14 -O3 -L./ -ldynamicboltz myprogram.cpp -o myprogram.o
```

## BMLA

To build the library called `bmla` - navigate to the dir `bmla`, then:
```
g++ -std=c++14 -O3 -c -fpic bmla.cpp ixn_param.cpp species.cpp
g++ -std=c++14 -O3 -shared -o libbmla.so bmla.o ../general.o species.o ixn_param.o ../lattice.o
```

Link against it using
```
g++ -std=c++14 -O3 -L./ -lbmla myprogram.cpp -o myprogram.o
```