
.. _program_listing_file_main_rossler_varying_ic.cpp:

Program Listing for File main_rossler_varying_ic.cpp
====================================================

- Return to documentation for :ref:`file_main_rossler_varying_ic.cpp`

.. code-block:: cpp

   #include "dynamic_boltzmann.hpp"
   #include "general.hpp"
   #include <iostream>
   
   using namespace DynamicBoltzmann;
   
   int main() {    
   
       // Dimensions
       std::vector<Dim> dims;
       dims.push_back(Dim("hA",DimType::H,"A",{"hA","hB","hC"},-2.0,0.5,20,-1.0));
       dims.push_back(Dim("hB",DimType::H,"B",{"hA","hB","hC"},-2.0,0.5,20,-1.0));
       dims.push_back(Dim("hC",DimType::H,"C",{"hA","hB","hC"},-2.0,0.5,20,-1.0));
       dims.push_back(Dim("jAA",DimType::J,"A","A",{"jAA","jAB","jAC"},-1.5,1.5,20,0.0));
       dims.push_back(Dim("jAB",DimType::J,"A","B",{"jAA","jAB","jBB"},-1.5,1.5,20,0.0));
       dims.push_back(Dim("jAC",DimType::J,"A","C",{"jAA","jAC","jCC"},-1.5,1.5,20,0.0));
       dims.push_back(Dim("jBB",DimType::J,"B","B",{"jAA","jAB","jBB"},-1.5,1.5,20,0.0));
       dims.push_back(Dim("jBC",DimType::J,"B","C",{"jBB","jBC","jCC"},-1.5,1.5,20,0.0));
       dims.push_back(Dim("jCC",DimType::J,"C","C",{"jAC","jBC","jCC"},-1.5,1.5,20,0.0));
   
       // Times
       double t_max=0.8;
       int n_t = 80;
   
       // Opt params
       int batch_size = 50;
       int n_annealing = 3000;
       int box_length = 10;
       double dopt = 0.1;
       int n_opt = 250;
       
       // Opt problem
       OptProblem opt(dims,{"A","B","C"},t_max,n_t,batch_size,n_annealing,box_length,dopt,n_opt);
   
       // Print out validation
       // opt.validate_setup();
   
       // Add filenames
       for (int i=0; i<100; i++) {
           opt.add_fname("rossler_varying_ic/lattice_v" + pad_str(i,2) + "/lattice/");
       };
   
       // Solve
       std::cout << "Solving..." << std::endl;
       opt.solve_varying_ic(true);
       std::cout << "fin." << std::endl;
   
       return 0;
   }
