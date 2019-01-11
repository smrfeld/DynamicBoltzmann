// vector
#ifndef VECTOR_h
#define VECTOR_h
#include <vector>
#endif

// list
#ifndef LIST_h
#define LIST_h
#include <list>
#endif

// string
#ifndef STRING_h
#define STRING_h
#include <string>
#endif

// map
#ifndef MAP_h
#define MAP_h
#include <map>
#endif

// Diagnostic flags
#define DIAG_DIFFUSE 0

/************************************
* Namespace for LatticeGillespie
************************************/

namespace LatticeGillespie {

	/****************************************
	Main simulation class
	****************************************/

	class Simulation
	{
	private:

		class Impl;
		std::unique_ptr<Impl> _impl;

	public:

		/********************
		Constructor/Destructor
		********************/

		Simulation(double dt, int box_length, int dim=3);
		Simulation(Simulation&& other); // movable but no copies
	    Simulation& operator=(Simulation&& other); // movable but no copies
		~Simulation();

		/********************
		Add species
		********************/

		void add_species(std::string name, bool conserved = false);

		/********************
		 Add a reaction
		********************/

		// Unimolecular rxn
		void add_uni_rxn(std::string name, double kr, std::string r);
		void add_uni_rxn(std::string name, double kr, std::string r, std::string p);
		void add_uni_rxn(std::string name, double kr, std::string r, std::string p1, std::string p2);

		// Bimolecular rxn
		void add_bi_rxn(std::string name, double prob, std::string r1, std::string r2);
		void add_bi_rxn(std::string name, double prob, std::string r1, std::string r2, std::string p);
		void add_bi_rxn(std::string name, double prob, std::string r1, std::string r2, std::string p1, std::string p2);

		/********************
		Populate lattice
		********************/

		void populate_lattice(std::map<std::string,int> counts);
		void populate_lattice(std::map<std::string,double> &h_dict, std::map<std::string,std::map<std::string,double>> &j_dict, int n_steps);
		void populate_lattice(std::map<std::string,double> &h_dict, std::map<std::string,std::map<std::string,double>> &j_dict, std::map<std::string, std::map<std::string,std::map<std::string,double>>> &k_dict, int n_steps);

		/********************
		Run simulation
		********************/

		void run(int n_timesteps, bool verbose = true, bool write_counts = false, bool write_nns = false, bool write_latt = false, int write_step = 20, int write_version_no = 0);

		/********************
		Write lattice
		********************/

		void write_lattice(int index, int write_version_no);
	};
};