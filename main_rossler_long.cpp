#include "dynamic_boltzmann.hpp"
#include "general.hpp"
#include <iostream>

using namespace DynamicBoltzmann;

int main() {	

	// Dimensions
	std::vector<Dim> dims;
	dims.push_back(Dim("hA",DimType::H,"A",{"hA","hB","hC"},-1.5,0.0,20,-1.0));
	dims.push_back(Dim("hB",DimType::H,"B",{"hA","hB","hC"},-1.1,0.3,20,-1.0));
	dims.push_back(Dim("hC",DimType::H,"C",{"hA","hB","hC"},-2.0,0.0,20,-1.0));
	dims.push_back(Dim("jAA",DimType::J,"A","A",{"jAA","jAB","jAC"},-0.1,1.4,20,0.0));
	dims.push_back(Dim("jAB",DimType::J,"A","B",{"jAA","jAB","jBB"},-0.3,0.8,20,0.0));
	dims.push_back(Dim("jAC",DimType::J,"A","C",{"jAA","jAC","jCC"},-0.8,0.3,20,0.0));
	dims.push_back(Dim("jBB",DimType::J,"B","B",{"jAA","jAB","jBB"},-0.1,1.1,20,0.0));
	dims.push_back(Dim("jBC",DimType::J,"B","C",{"jBB","jBC","jCC"},-0.5,0.8,20,0.0));
	dims.push_back(Dim("jCC",DimType::J,"C","C",{"jAC","jBC","jCC"},-0.25,1.5,20,0.0));

	// Opt params 
	int batch_size = 5;
	int n_annealing = 100000;
	int box_length = 10;
	double dopt = 0.01;
	int n_opt = 200;
	
	double t_max;
	int n_t;
	int n_t_old;

	/********************
	Solve all at once
	********************/

	t_max = 2.0;
	n_t = 200;

	// Opt problem
	OptProblem *opt1 = new OptProblem(dims,{"A","B","C"},t_max,n_t,batch_size,n_annealing,box_length,dopt,n_opt);

	// Write only last bf
	opt1->set_flag_write_bf_only_final();

	// Add filenames
	for (int i=0; i<300; i++) {
		opt1->add_fname("rossler_long/lattice_v" + pad_str(i,2) + "/lattice/");
	};

	// Set dir
	opt1->set_dir_io("rossler_long_data/");

	// Solve
	std::cout << "Solving..." << std::endl;
	opt1->solve(true);
	std::cout << "fin." << std::endl;

	/********************
	Solve in segments
	********************/

	/*
	// Times
	t_max=0.0;
	n_t = 0;

	for (int i_loop=1; i_loop<=4; i_loop++)
	{
		n_t_old = n_t;
		t_max += 0.2;
		n_t += 20;

		if (i_loop==1) {
			dopt = 0.1;
		} else {
			dopt = 0.025;
		};

		// Opt problem
		OptProblem *opt1 = new OptProblem(dims,{"A","B","C"},t_max,n_t,batch_size,n_annealing,box_length,dopt,n_opt);

		// Write only last bf
		opt1->set_flag_write_bf_only_final();

		// Add filenames
		for (int i=0; i<1; i++) {
			opt1->add_fname("rossler_long/lattice_v" + pad_str(i,2) + "/lattice/");
		};

		// Set dir
		opt1->set_dir_io("rossler_long_data_"+pad_str(i_loop,1)+"/");

		// Read in the basis funcs
		if (i_loop!=0) {
			opt1->read_bf("F_hA","rossler_long_data_"+pad_str(i_loop-1,1)+"/F/F_hA_"+pad_str(n_t_old,4)+".txt");	
			opt1->read_bf("F_hB","rossler_long_data_"+pad_str(i_loop-1,1)+"/F/F_hB_"+pad_str(n_t_old,4)+".txt");	
			opt1->read_bf("F_hC","rossler_long_data_"+pad_str(i_loop-1,1)+"/F/F_hC_"+pad_str(n_t_old,4)+".txt");	
			opt1->read_bf("F_jAA","rossler_long_data_"+pad_str(i_loop-1,1)+"/F/F_jAA_"+pad_str(n_t_old,4)+".txt");	
			opt1->read_bf("F_jAB","rossler_long_data_"+pad_str(i_loop-1,1)+"/F/F_jAB_"+pad_str(n_t_old,4)+".txt");	
			opt1->read_bf("F_jAC","rossler_long_data_"+pad_str(i_loop-1,1)+"/F/F_jAC_"+pad_str(n_t_old,4)+".txt");	
			opt1->read_bf("F_jBB","rossler_long_data_"+pad_str(i_loop-1,1)+"/F/F_jBB_"+pad_str(n_t_old,4)+".txt");	
			opt1->read_bf("F_jBC","rossler_long_data_"+pad_str(i_loop-1,1)+"/F/F_jBC_"+pad_str(n_t_old,4)+".txt");	
			opt1->read_bf("F_jCC","rossler_long_data_"+pad_str(i_loop-1,1)+"/F/F_jCC_"+pad_str(n_t_old,4)+".txt");	
		};

		// Solve
		std::cout << "Solving..." << std::endl;
		opt1->solve(true);
		std::cout << "fin." << std::endl;

		// Clean up
		delete opt1;

	};
	*/
	
	return 0;
}