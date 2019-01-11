
.. _program_listing_file_main_annihilation.cpp:

Program Listing for File main_annihilation.cpp
==============================================

- Return to documentation for :ref:`file_main_annihilation.cpp`

.. code-block:: cpp

   #include "dynamic_boltzmann.hpp"
   #include "general.hpp"
   #include <iostream>
   
   using namespace DynamicBoltzmann;
   
   int main() {
   
       /********************
       Test optimization problem
       ********************/
   
       // Dimensions vec
       std::vector<Dim> dims;
       dims.push_back(Dim("h",DimType::H,"A",{"h","J"},-0.6,0.55,21,0.5));
       dims.push_back(Dim("J",DimType::J,"A","A",{"h","J"},-0.7,0.15,21,0.1));
   
       // Times
       double t_max=1.0;
       int n_t = 101;
   
       // Opt params
       int batch_size = 10;
       int n_annealing = 500;
       int box_length = 10;
       double dopt = 0.005;
       int n_opt = 100;
       
       // Init
       OptProblem opt(dims, {"A"}, t_max, n_t, batch_size, n_annealing, box_length, dopt, n_opt);
   
       // Add filenames
       for (int i=0; i<100; i++) {
           opt.add_fname("annihilation/lattice_v" + pad_str(i,2) + "/lattice/");
       };
   
       // Validate setup
       opt.validate_setup();
   
       // Solve
       std::cout << "Solving..." << std::endl;
       opt.solve(true);
       std::cout << "fin." << std::endl;
   
       return 0;
   };
