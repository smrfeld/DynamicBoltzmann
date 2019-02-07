# Dynamic Boltzmann Machine for Lattice Chemical Kinetics

## Dependencies

`armadillo` library.

`q3c1` library for the Q3 C1 finite elements. Get it from [here](https://github.com/smrfeld/Q3-C1-Finite-Elements).

## Installation

Use `cmake`:
```
mkdir build
cd build
cmake ..
make
make install
```
The default install locations are `/usr/local/include` and `/usr/local/lib`.

If not installing, make sure to move the library to somewhere it will be found, or at worst
```
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/absolute/path/to/lib
```
In this case, you will also need to specify using `-L` the location of the library at linking time (`-L` is needed for finding the library at linking time; `DYLD_LIBRARY_PATH` is at runtime).

## Usage

To access the headers, simply use: `include <dblz>` or `include <bmla>. 

Link your program using:
```
g++ -std=c++14 -O3 -ldblz -lq3c1 myprogram.cpp -o myprogram.o
g++ -std=c++14 -O3 -lbmla myprogram.cpp -o myprogram.o
```
(this assumes you have install `-lq3c1`).

## Namespace

The namespace is `dblz` ("DynamicBoltzmann").

## Technical notes (POSSIBLY OUTDATED!)

### Binary vs probabilistic units

The guidelines for RBMs in Geoffrey Hintons guide "A Practical Guide to Training Restricted Boltzmann Machines" state:
* "When the hidden units are being driven by data, always use stochastic binary states. When they are being driven by reconstructions, always use probabilities without sampling."

* "Assuming the visible units use the logistic function, use real-valued probabilities for both the data and the reconstructions."

* "When collecting the pairwise statistics for learning weights or the individual statistics for learning biases, use the probabilities, not the binary states, and make sure the weights have random initial values to break symmetry."

These rules lead to problems here - particularly the first rule: **"When they are being driven by reconstructions, always use probabilities without sampling."** Consider the following:

1. The stochastic simulations are binary by nature. 
2. Assume we activate the hidden layer from these by also using binary units. Hence, when evaluating the **awake phase** moment `<v * h>` between visible and hidden units will have many terms of the form `0*1` or `1*0` or `0*0`, leading to many null contributions to `sum_i <v_i * h_i>`.
3. Now we switch do using probabilistic units for doing the sampling for the **asleep phase** moment. In this case, in the contributions from before, we replace `0` by `some small #`. The upshot is that where contributions were null before, here they are **greater than zero.**

Therefore, the **asleep moment will always be larger than the awake phase moment.** This drives the weights `W` to `-infinity`!

Instead, we simply use **binary units everywhere**, i.e.:
1. The stochastic simulations are binary by nature. 
2. Use binary units for the hidden units for the first round of sampling (driven by data). We then compute the awake phase moment.
3. Use binary units for other rounds of sampling both hidden and visible units, and to evaluate the asleep phase moments.

This avoids the `infinity` problem.

# BMLA - To learn initial conditions

To learn initial conditions, see the `bmla` folder [here](bmla).