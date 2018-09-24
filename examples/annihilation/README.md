# Simple unimolecular annihilation example A->0

## To generate data

NOTE: To be update soon!

To generate data via stochastic simulations:

In the `stoch_sim` directory:

Run `python create_dirs.py` to create the directory structures.

Modify `main.py` as needed. Then compile using `g++ -std=c++14 -L./ -lgillespie3d main.cpp -o main.o` and run it.