#ifndef LATTICE_h
#define LATTICE_h
#include "lattice.hpp"
#endif

#ifndef SPECIES_H
#define SPECIES_H
#include "species.hpp"
#endif

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	class BMLA {

	private:

		// The lattice to learn
		Lattice _latt;

		// Number of annealing steps
		int _n_annealing;

		// The species present
		std::list<Species> _species;

		// The initial guesses at the couplings
		std::vector<Species*> _h_species; 
		std::vector<std::pair<Species*,Species*>> _j_species; 
		std::vector<double> _h0;
		std::vector<double> _j0;

		// The current h,j
		// For reasons unknown (poor code planning), first dimension just has one entry
		double **_soln_traj;
		int _t_fake;
		int _n_params;

		// Initial and final moments
		double *_moms_init;
		double *_moms_final;

		// Optimization step size
		double _dopt;

		// No opt steps
		int _n_opt;

		// Batch size
		int _n_batch;

	public:

		// Constructor
		// Note: j0 should not be doubly linked
		BMLA(int box_length, int n_annealing, double dopt, std::vector<std::string> species, std::map<std::string,double> h0, std::map<std::string,std::map<std::string,double>> j0, int n_params, int n_opt, int n_batch);
		~BMLA();

		// Solve for the h,j corresponding to a given lattice
		void solve(std::string fname, bool verbose=false);

		// Update the initial params
		void update_initial_params(std::map<std::string,double> h0, std::map<std::string,std::map<std::string,double>> j0);

		// Write out the solutions
		void write(std::string fname);

	};

};


