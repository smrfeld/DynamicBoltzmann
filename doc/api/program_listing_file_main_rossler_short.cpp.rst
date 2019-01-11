
.. _program_listing_file_main_rossler_short.cpp:

Program Listing for File main_rossler_short.cpp
===============================================

- Return to documentation for :ref:`file_main_rossler_short.cpp`

.. code-block:: cpp

   #include "dynamic_boltzmann.hpp"
   #include "general.hpp"
   #include <iostream>
   
   using namespace DynamicBoltzmann;
   
   int main() {    
   
       // Dimensions
       std::vector<Dim> dims;
       dims.push_back(Dim("hA",DimType::H,"A",{"hA","hB","hC"},-1.6,-0.3,20,-1.0));
       dims.push_back(Dim("hB",DimType::H,"B",{"hA","hB","hC"},-1.1,0.8,20,-1.0));
       dims.push_back(Dim("hC",DimType::H,"C",{"hA","hB","hC"},-1.1,-0.7,20,-1.0));
       dims.push_back(Dim("jAA",DimType::J,"A","A",{"jAA","jAB","jAC"},-0.1,1.4,20,0.0));
       dims.push_back(Dim("jAB",DimType::J,"A","B",{"jAA","jAB","jBB"},-0.3,0.5,20,0.0));
       dims.push_back(Dim("jAC",DimType::J,"A","C",{"jAA","jAC","jCC"},-0.4,0.3,20,0.0));
       dims.push_back(Dim("jBB",DimType::J,"B","B",{"jAA","jAB","jBB"},-0.1,1.1,20,0.0));
       dims.push_back(Dim("jBC",DimType::J,"B","C",{"jBB","jBC","jCC"},-0.1,0.4,20,0.0));
       dims.push_back(Dim("jCC",DimType::J,"C","C",{"jAC","jBC","jCC"},-0.1,1.5,20,0.0));
   
       // Times
       double t_max=0.4;
       int n_t = 40;
   
       // Opt params
       int batch_size = 5;
       int n_annealing = 3000;
       int box_length = 10;
       double dopt = 0.1;
       int n_opt = 250;
       
       // Opt problem
       OptProblem opt(dims,{"A","B","C"},t_max,n_t,batch_size,n_annealing,box_length,dopt,n_opt);
   
       // Print out validation
       // opt.validate_setup();
   
       // Add filenames
       for (int i=0; i<1; i++) {
           opt.add_fname("rossler_long/lattice_v" + pad_str(i,2) + "/lattice/");
       };
   
       // Solve
       std::cout << "Solving..." << std::endl;
       opt.solve(true);
       std::cout << "fin." << std::endl;
   
       return 0;
   }
