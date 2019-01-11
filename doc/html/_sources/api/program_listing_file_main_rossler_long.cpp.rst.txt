
.. _program_listing_file_main_rossler_long.cpp:

Program Listing for File main_rossler_long.cpp
==============================================

- Return to documentation for :ref:`file_main_rossler_long.cpp`

.. code-block:: cpp

   #include "dynamic_boltzmann.hpp"
   #include "general.hpp"
   #include <iostream>
   
   using namespace DynamicBoltzmann;
   
   int main() {    
   
       // Dimensions
       std::vector<Dim> dims;
       dims.push_back(Dim("hA",DimType::H,"A",{"hA","hB","hC"},-2.0,-0.3,20,-1.0));
       dims.push_back(Dim("hB",DimType::H,"B",{"hA","hB","hC"},-1.1,0.8,20,-1.0));
       dims.push_back(Dim("hC",DimType::H,"C",{"hA","hB","hC"},-1.1,0.2,20,-1.0));
       dims.push_back(Dim("jAA",DimType::J,"A","A",{"jAA","jAB","jAC"},-0.1,1.8,20,0.0));
       dims.push_back(Dim("jAB",DimType::J,"A","B",{"jAA","jAB","jBB"},-0.3,0.5,20,0.0));
       dims.push_back(Dim("jAC",DimType::J,"A","C",{"jAA","jAC","jCC"},-0.8,0.3,20,0.0));
       dims.push_back(Dim("jBB",DimType::J,"B","B",{"jAA","jAB","jBB"},-0.5,1.1,20,0.0));
       dims.push_back(Dim("jBC",DimType::J,"B","C",{"jBB","jBC","jCC"},-0.4,0.4,20,0.0));
       dims.push_back(Dim("jCC",DimType::J,"C","C",{"jAC","jBC","jCC"},-0.1,1.5,20,0.0));
   
       // Opt params 
       int batch_size = 5;
       int n_annealing = 3000;
       int box_length = 10;
       double dopt = 0.2;
       int n_opt = 200;
       
       double t_max;
       int n_t;
   
       /********************
       Solve part 1
       ********************/
   
       // Times
       t_max=0.4;
       n_t = 40;
   
       // Opt problem
       OptProblem *opt1 = new OptProblem(dims,{"A","B","C"},t_max,n_t,batch_size,n_annealing,box_length,dopt,n_opt);
   
       // Add filenames
       for (int i=0; i<1; i++) {
           opt1->add_fname("rossler_long/lattice_v" + pad_str(i,2) + "/lattice/");
       };
   
       // Set dir
       opt1->set_dir_io("rossler_long_data_1/");
   
       // Solve
       std::cout << "Solving..." << std::endl;
       opt1->solve(true);
       std::cout << "fin." << std::endl;
   
       // Clean up
       delete opt1;
   
       /********************
       Solve part 2
       ********************/
   
       // Times
       t_max=0.8;
       n_t = 80;
   
       // Opt problem
       OptProblem *opt2 = new OptProblem(dims,{"A","B","C"},t_max,n_t,batch_size,n_annealing,box_length,dopt,n_opt);
   
       // Add filenames
       for (int i=0; i<1; i++) {
           opt2->add_fname("rossler_long/lattice_v" + pad_str(i,2) + "/lattice/");
       };
   
       // Set dir
       opt2->set_dir_io("rossler_long_data_2/");
   
       // Read in the basis funcs
       opt2->read_bf("F_hA","rossler_long_data_1/F/F_hA_0200.txt");    
       opt2->read_bf("F_hB","rossler_long_data_1/F/F_hB_0200.txt");    
       opt2->read_bf("F_hC","rossler_long_data_1/F/F_hC_0200.txt");    
       opt2->read_bf("F_jAA","rossler_long_data_1/F/F_jAA_0200.txt");  
       opt2->read_bf("F_jAB","rossler_long_data_1/F/F_jAB_0200.txt");  
       opt2->read_bf("F_jAC","rossler_long_data_1/F/F_jAC_0200.txt");  
       opt2->read_bf("F_jBB","rossler_long_data_1/F/F_jBB_0200.txt");  
       opt2->read_bf("F_jBC","rossler_long_data_1/F/F_jBC_0200.txt");  
       opt2->read_bf("F_jCC","rossler_long_data_1/F/F_jCC_0200.txt");  
   
       // Solve
       std::cout << "Solving..." << std::endl;
       opt2->solve(true);
       std::cout << "fin." << std::endl;
   
       return 0;
   }
