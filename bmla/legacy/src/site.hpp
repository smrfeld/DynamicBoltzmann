// Species
#ifndef SPECIES_h
#define SPECIES_h
#include "species.hpp"
#endif

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	General functions
	****************************************/

	class Site;
	typedef std::list<Site> lattice;
	typedef std::list<Site>::iterator latt_it;

	/****************************************
	Structure to hold iterators to neighbors
	****************************************/

	struct LattIt2 {
		latt_it lit1,lit2;
		LattIt2(latt_it l1, latt_it l2) {
			lit1 = l1;
			lit2 = l2;
		};
	};

	struct LattIt3 {
		latt_it lit1,lit2,lit3;
		LattIt3(latt_it l1, latt_it l2, latt_it l3) {
			lit1 = l1;
			lit2 = l2;
			lit3 = l3;
		};
	};

	/****************************************
	Class to hold a lattice site
	****************************************/

	class Site {
	private:

		// Probs of species/empty
		std::map<Species*, double> _probs;
		double _prob_empty;

		// Add/Remove counts on a given species (not _sp_binary, unless it is passed)
		void _remove_counts_on_species(Species *sp, double prob);
		void _add_counts_on_species(Species *sp, double prob);

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const Site& other);

	public:

		int dim;
		int x;
		int y;
		int z;
		std::vector<latt_it> nbrs;
		std::vector<LattIt2> nbrs_triplets;
		std::vector<LattIt3> nbrs_quartics;

		// Connectivity to any hidden units
		// A species-dependent graph :)
		std::vector<ConnectionVH*> hidden_conns;

		// Constructor
		Site(int xIn);
		Site(int xIn, int yIn);
		Site(int xIn, int yIn, int zIn);
		Site(const Site& other);
		Site(Site&& other);
		Site& operator=(const Site& other);
		Site& operator=(Site&& other);
		~Site();

		// Add a species possibility
		void add_species_possibility(Species* sp);

		// Get a probability
		// If binary, returns 1.0 for a certain species
		// Pass nullptr to get prob of empty
		double get_prob(Species *sp) const;
		// Get all probs
		const std::map<Species*, double>& get_probs() const;
		// Set probability
		// Pass nullptr to set probability of being empty
		void set_prob(Species *sp, double prob);
		// Set site to be empty
		void set_site_empty();
		// Set site to have binary probs
		void set_site_binary(Species *sp);

		// Get J and K activations
		// Go through all possible species: (probabilities) times (J of the coupling for the given species)
		double get_act_j(Species *sp) const;
		double get_act_k(Species *sp) const;

		// Check if site is empty
		bool empty() const;

		// Binarize the site
		void binarize();
	};
	// Comparator
	bool operator <(const Site& a, const Site& b);
	bool operator==(const Site& a, const Site& b);
	std::ostream& operator<<(std::ostream& os, const Site& s);

};