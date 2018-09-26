# Bimolecular Annihilation Example: A+B->0

To generate data, see `stoch_sim` folder.

To build, run `g++ -std=c++14 -ldynamicboltz -ldcubic main.cpp -o main.o`.

Before running, make the following dirs: `data_learned, data_learned/diff_eq_rhs, data_learned/moments, data_learned/adjoint, data_learned/ixn_params`.

Then run `./main.o`.

Analyze using the provided Mathematica files.