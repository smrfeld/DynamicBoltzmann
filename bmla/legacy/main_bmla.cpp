#include <iostream>
#include "bmla.hpp"
#include "general.hpp"
#include <fstream>

using namespace DynamicBoltzmann;

void read(std::string fname, std::map<std::string,double> &h0, std::map<std::string,std::map<std::string,double>> &j0) {

	std::ifstream f;
	f.open(fname);
	char frag[100]; // fragments of the line
	int i_frag=0;
	int i_line=0;
	std::string s1="",s2="";
	std::string val="";

	if (f.is_open()) { // make sure we found it
		while (!f.eof()) {
		    f >> frag;
		    if (i_line < 3) {
		    	if (i_frag==0) { 
		    		s1 += frag; i_frag++;
		    	} else if (i_frag==1) {
		    		val += frag;
		    		// Add to map
		    		h0[s1] = std::stod(val);
		    		// Reset
		    		i_line++;
		    		i_frag=0;
		    		s1="";s2="";val="";
		    	};
		    } else {
		    	if (i_frag==0) { 
		    		s1 += frag; i_frag++;
		    	} else if (i_frag==1) { 
		    		s2 += frag; i_frag++;
		    	} else if (i_frag==2) {
		    		val += frag;
		    		// Add to map
		    		j0[s1][s2] = std::stod(val);
		    		// Reset
		    		i_line++;
		    		i_frag=0;
		    		s1="";s2="";val="";
		    	};
		    };
		};
	};

	f.close();
};

int main() {

	int box_length = 10;
	int n_annealing = 1000000;
	double dopt = 0.001;

	std::vector<std::string> species = {"A","B","C"};

	std::map<std::string,double> h0;
	std::map<std::string,std::map<std::string,double>> j0;
	// Guess
	h0["A"] = -0.406;
	h0["B"] = -0.647;
	h0["C"] = -0.507;
	j0["A"]["A"] = 0.05;
	j0["A"]["B"] = -0.235;
	j0["A"]["C"] = -0.088;
	j0["B"]["B"] = 0.246;
	j0["B"]["C"] = -0.131;
	j0["C"]["C"] = -0.262;

	int n_params = 9;

	int n_opt = 10;

	int n_batch = 10;

	BMLA bmla(box_length,n_annealing,dopt,species,h0,j0,n_params,n_opt,n_batch);

	// Iterate over the things
	for (int i=0; i<142; i++) {
		std::cout << "--- File: " << i << " / 141 ---" << std::endl; 

		// Read in previous initial params
		read("../rossler_varying_ic/"+pad_str(i,3)+"/init.txt",h0,j0);

		// Update
		bmla.update_initial_params(h0,j0);

		// Solve
		bmla.solve("../rossler_varying_ic/"+pad_str(i,3)+"/0000.txt",true);
		
		// Write
		bmla.write("../rossler_varying_ic/"+pad_str(i,3)+"/init.txt");
	};

	return 0;
};