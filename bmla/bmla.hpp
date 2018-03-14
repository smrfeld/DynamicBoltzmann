#ifndef LATTICE_h
#define LATTICE_h
#include "../lattice.hpp"
#endif

#ifndef IXN_PARAM_h
#define IXN_PARAM_h
#include "ixn_param.hpp"
#endif

#ifndef LIST_H
#define LIST_H
#include <list>
#endif

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Struct for specifying a dimension
	****************************************/

	// Type of dimension
	enum DimType { H, J };

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

		// Number of dimensions
		int _n_param;

		// List of interaction parameters
		std::list<IxnParam> _ixn_params;

		// Species present
		std::list<Species> _species;

		// Batch size
		int _n_batch;

		// Number of annealing steps
		int _n_annealing;

		// The lattice to learn
		Lattice _latt;

		// Optimization step size
		double _dopt;

		// No opt steps
		int _n_opt;

		// If the MSE dips below this, quit
		bool _mse_quit_mode;
		double _mse_quit; // in percent

		// L2 reg
		bool _l2_reg;
		double _lambda;

		// Print
		void _print_ixn_params(bool new_line=true) const;
		void _print_moments() const;
		void _print_mse(bool new_line=true) const;

		// Get the mse
		double _get_mse() const;

		// Search functions
		Species* _find_species(std::string name);
		IxnParam* _find_ixn_param(std::string name);

		// Constructor helpers
		void _clean_up();
		void _copy(const BMLA& other);
		void _copy(BMLA&& other);

	public:

		// Constructor
		BMLA(std::vector<Dim> dims, std::vector<std::string> species, int batch_size, int n_annealing, int box_length, double dopt, int n_opt, int lattice_dim=3);
		BMLA(const BMLA& other);
		BMLA(BMLA&& other);
		BMLA& operator=(const BMLA& other);
	    BMLA& operator=(BMLA&& other);
		~BMLA();

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


