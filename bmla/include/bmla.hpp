#ifndef LIST_H
#define LIST_H
#include <list>
#endif

#ifndef STRING_H
#define STRING_H
#include <string>
#endif

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Struct for specifying a dimension
	****************************************/

	// Type of dimension
	enum DimType { H, J, W };

	struct Dim {
		// Name
		std::string name;

		// Type
		DimType type;

		// Name of associated species
		std::string species1;
		std::string species2;

		// Guess
		double guess;

		// Constructor
		Dim(std::string name, DimType type, std::string species, double guess);
		Dim(std::string name, DimType type, std::string species1, std::string species2, double guess);
	};

	/****************************************
	BMLA
	****************************************/

	class BMLA {

	private:

		class Impl;
		std::unique_ptr<Impl> _impl;

	public:

		// Constructor
		BMLA(std::vector<Dim> dims, std::vector<std::string> species, int batch_size, int box_length, double dopt, int n_opt, int lattice_dim=3);
		BMLA(BMLA&& other); // movable but no copies
	    BMLA& operator=(BMLA&& other); // movable but no copies
		~BMLA();

		// Any dim
		void add_hidden_unit(std::vector<std::vector<int>> lattice_idxs, std::string species);
		// 1D specific
		void add_hidden_unit(std::vector<int> lattice_idxs, std::string species);

		// Set the number of CD steps (default = 1)
		void set_n_cd_steps(int n_steps);

		// Set and turn on l2 regularizer
		void set_l2_reg(double lambda);

		// Set and turn on MSE quit mode
		void set_mse_quit(double mse_quit);

		// Solve for the h,j corresponding to a given lattice
		void solve(std::string fname, bool verbose=false);

		// Update the initial params
		void read(std::string fname);

		// Write out the solutions
		void write(std::string fname);
	};
};


