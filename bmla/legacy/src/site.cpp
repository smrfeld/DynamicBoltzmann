#include "site.hpp"
#include <iostream>
#include <fstream>
#include <numeric>
#include "../include/general.hpp"
#include "math.h"
#include <ctime>
#include <sstream>
#include <random>

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Struct to hold a lattice Site
	****************************************/

	// Constructor
	Site::Site(int xIn) : Site(xIn,0,0) { dim=1; };
	Site::Site(int xIn, int yIn) : Site(xIn,yIn,0) { dim=2; };
	Site::Site(int xIn, int yIn, int zIn) {
		dim=3;
		x = xIn;
		y = yIn;
		z = zIn;
		_prob_empty = 1.0; // default = empty
	};	
	Site::Site(const Site& other) {
		_copy(other);
	};
	Site::Site(Site&& other) {
		_copy(other);
		other._reset();
	};
	Site& Site::operator=(const Site& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	Site& Site::operator=(Site&& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
	};
	Site::~Site() {
		_clean_up();
	};
	void Site::_clean_up() {
		// Nothing....
	};
	void Site::_reset() {
		dim = 0;
		x = 0;
		y = 0;
		z = 0;
		nbrs.clear();
		nbrs_triplets.clear();
		nbrs_quartics.clear();
		hidden_conns.clear();
		_prob_empty = 0.0;
		_probs.clear();
	};
	void Site::_copy(const Site& other) {
		dim = other.dim;
		x = other.x;
		y = other.y;
		z = other.z;
		nbrs = other.nbrs;
		nbrs_triplets = other.nbrs_triplets;
		nbrs_quartics = other.nbrs_quartics;
		hidden_conns = other.hidden_conns;
		_prob_empty = other._prob_empty;
		_probs = other._probs;
	};

	// Comparator
	bool operator <(const Site& a, const Site& b) {
		if (a.dim == 1) {
	    	return a.x < b.x;
	    } else if (a.dim == 2) {
	    	return std::tie(a.x, a.y) < std::tie(b.x, b.y);
	    } else if (a.dim == 3) {
	    	return std::tie(a.x, a.y, a.z) < std::tie(b.x, b.y, b.z);
	    } else {
	    	return false;
	    };
	};
	bool operator==(const Site& a, const Site& b) {
		if (a.dim == 1) {
	    	return a.x == b.x;
	    } else if (a.dim == 2) {
	    	return std::tie(a.x, a.y) == std::tie(b.x, b.y);
	    } else if (a.dim == 3) {
	    	return std::tie(a.x, a.y, a.z) == std::tie(b.x, b.y, b.z);
	    } else {
	    	return false;
	    };
	}; 
	std::ostream& operator<<(std::ostream& os, const Site& s)
	{
		if (s.empty()) {
			return os;
		};

		if (s.dim == 1) {
		    os << s.x << " ";
		} else if (s.dim == 2) {
		    os << s.x << " " << s.y << " ";
		} else if (s.dim == 3) {
		    os << s.x << " " << s.y << " " << s.z << " ";
		};

		// Get probs
		const std::map<Species*, double> prs = s.get_probs();

		// All probs
		for (auto pr: prs) {
			if (pr.second > 0.0) {
				os << pr.first->name() << " " << pr.second << " ";	
			};
		};

	    return os;
    };

	// Add a species possibility
	void Site::add_species_possibility(Species* sp) {
		_probs[sp] = 0.0;
	};

	// Get a probability
	// nullptr for empty
	double Site::get_prob(Species *sp) const {
		// nullptr for prob of empty
		if (sp == nullptr) {
			return _prob_empty;
		};

		auto it = _probs.find(sp);
		if (it != _probs.end()) {
			return it->second;
		} else {
			return 0.0;
		};
	};

	// Get all probs (excluding prob of empty)
	const std::map<Species*, double>& Site::get_probs() const {
		return _probs;
	};

	// Set probability
	// Pass nullptr to set probability of being empty
	void Site::set_prob(Species *sp, double prob) {

		// nullptr for prob of empty
		if (sp == nullptr) {
			_prob_empty = prob;
			return;
		};

		// Remove the old counts
		if (_probs[sp] > 0.0) {
			_remove_counts_on_species(sp,_probs[sp]);
		};

		// Store
		// std::cout << "Set prob for site: " << x << " species " << sp->name() << " to: " << prob << std::endl;
		_probs[sp] = prob;

		// Increment counts on species
		_add_counts_on_species(sp,prob);
	};

	// Set a site to be empty
	void Site::set_site_empty() {
		// Remove counts on existing species
		for (auto pr: _probs) {
			if (pr.second > 0.0) {
				_remove_counts_on_species(pr.first,pr.second);
				// Clear
				_probs[pr.first] = 0.0;
			};
		};

		// Empty prob = 1
		_prob_empty = 1.0;
	};

	// Set site to have binary probs
	void Site::set_site_binary(Species *sp) {
		// Clear
		set_site_empty();
		// Make
		set_prob(sp,1.0);
	};

	// Get J and K activations
	// Go through all possible species probabilities x J of the coupling for the given species
	double Site::get_act_j(Species *sp) const {
		double act=0.0;
		// Go through all nbrs
		for (auto lit: nbrs) {
			// Get all probs
			const std::map<Species*, double> prs = lit->get_probs();
			// Go through all probs
			for (auto pr: prs) {
				// J * prob
				act += sp->j(pr.first) * pr.second;
			};
		};
		return act;
	};
	double Site::get_act_k(Species *sp) const {
		double act=0.0;
		// Go through all pairs to consider
		for (auto trip: nbrs_triplets) {
			// Get all probs
			const std::map<Species*, double> prs1 = trip.lit1->get_probs();
			const std::map<Species*, double> prs2 = trip.lit2->get_probs();
			// Go through all probs
			for (auto pr1: prs1) {
				for (auto pr2: prs2) {
					// K * prob * prob
					act += sp->k(pr1.first,pr2.first) * pr1.second * pr2.second;
				};
			};
		};
		return act;
	};

	// Is site empty
	bool Site::empty() const {
		if (_prob_empty == 1.0) {
			return true;
		} else {
			return false;
		};
	};

	// Binarize the site
	void Site::binarize() {
		// Propensity vector
		std::vector<double> props;
		std::vector<Species*> sp_vec;
		props.push_back(0.0);
		props.push_back(_prob_empty);
		for (auto pr: _probs) {
			props.push_back(props.back()+pr.second);
			sp_vec.push_back(pr.first);
		};

		int i = sample_prop_vec(props);
		if (i==0) {
			set_site_empty();
		} else {
			set_site_binary(sp_vec[i-1]);
		};
	};

	/****************************************
	Site - PRIVATE
	****************************************/

	void Site::_remove_counts_on_species(Species *sp, double prob) {
		_add_counts_on_species(sp,-1.0*prob);
	};
	void Site::_add_counts_on_species(Species *sp, double prob) {

		// Counts
		sp->count_increment(prob);

		// NNs, if needed
		for (auto nbr_it: nbrs) {
			// Get all the probs
			const std::map<Species*, double> prs = nbr_it->get_probs();
			for (auto pr: prs) {
				// Increment
				sp->nn_count_increment(pr.first,prob*pr.second);
			};
		};

		// Triplets, if needed
		for (auto trip: nbrs_triplets) {
			// Get all the probs
			const std::map<Species*, double> prs1 = trip.lit1->get_probs();
			const std::map<Species*, double> prs2 = trip.lit2->get_probs();
			for (auto pr1: prs1) {
				for (auto pr2: prs2) {
					// Increment
					sp->triplet_count_increment(pr1.first,pr2.first,prob*pr1.second*pr2.second);
				};
			};
		};

		// Quartics, if needed
		for (auto quart: nbrs_quartics) {
			// Get all the probs
			const std::map<Species*, double> prs1 = quart.lit1->get_probs();
			const std::map<Species*, double> prs2 = quart.lit2->get_probs();
			const std::map<Species*, double> prs3 = quart.lit3->get_probs();
			for (auto pr1: prs1) {
				for (auto pr2: prs2) {
					for (auto pr3: prs3) {
						// Increment
						sp->quartic_count_increment(pr1.first,pr2.first,pr3.first,prob*pr1.second*pr2.second*pr3.second);
					};
				};
			};
		};
	};

};