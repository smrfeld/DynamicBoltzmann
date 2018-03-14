#ifndef VECTOR_h
#define VECTOR_h
#include <vector>
#endif

#ifndef STRING_h
#define STRING_h
#include <string>
#endif

#ifndef SPECIES_h
#define SPECIES_h
#include "species.hpp"
#endif

#ifndef GRID_h
#define GRID_h
#include "grid.hpp"
#endif

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Forward declare
	****************************************/

	class BasisFunc;

	/****************************************
	Interaction parameter
	****************************************/

	// Enumeration of type of dimension
	enum IxnParamType { Hp, Jp };

	class IxnParamTraj : public Grid {

	private:

		// Type
		IxnParamType _type;

		// Species
		Species *_sp1;
		Species *_sp2;

		// Number of time points in these trajs
		int _n_t;

		// Values over time
		double *_vals;

		// Initial value
		double _val0;

		// Awake and asleep moments over time
		double *_asleep;
		double *_awake;

		// Pointer to the basis function
		BasisFunc *_bf;

		// Copy, clean up
		void _copy(const IxnParamTraj& other);
		void _clean_up();

	public:

		// Constructor
		IxnParamTraj(std::string name, IxnParamType type, Species *sp, double min, double max, int n, double val0, int n_t);
		IxnParamTraj(std::string name, IxnParamType type, Species *sp1, Species *sp2, double min, double max, int n, double val0, int n_t);
		IxnParamTraj(const IxnParamTraj& other);
		IxnParamTraj & operator=(const IxnParamTraj& other);
		~IxnParamTraj();

		// Set/check basis func pointer
		void set_basis_func_ptr(BasisFunc* bf);
		bool is_bf(BasisFunc *bf);

		// Set IC
		void set_init_cond(double val);

		// Validate setup
		void validate_setup() const;

		// Getters/setters
		double get_at_time(int it) const;

		// Calculate the next step
		// Returns false if new point is out of grid
		bool calculate_at_time(int it_next, double dt);

		// Moments from lattice
		enum MomentType {AWAKE, ASLEEP};
		void moments_reset();
		void moments_retrieve_at_time(MomentType moment_type, int it);
		void moments_retrieve_at_time(MomentType moment_type, int it, int batch_size);
		double moments_diff_at_time(int it);

		// Write into an ofstream
		void write_vals(std::string dir, int idx, int n_t_traj) const;
		void write_vals(std::string dir, int idx1, int idx2, int n_t_traj) const;
		void write_moments(std::string dir, int idx, int n_t_traj) const;
		void write_moments(std::string dir, int idx1, int idx2, int n_t_traj) const;
	};
};

