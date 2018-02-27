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

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Dimension information
	****************************************/

	// Enumeration of type of dimension
	enum DimType { NONE, H, J };

	struct Dim
	{
		// Name
		std::string name;

		// Type = NONE, H, J
		DimType type;

		// Size = no pts
		int n;
		// Min/max/increment
		double min;
		double max;
		double delta;
		// Grid of vals
		double *grid;

		// Species involved
		std::string species1;
		std::string species2;
		Species *sp1;
		Species *sp2;

		// Awake/asleep moment
		double awake;
		double asleep;

		// Constructors
		Dim(std::string nameIn, double minIn, double maxIn, int nIn);
		Dim(std::string nameIn, double minIn, double maxIn, int nIn, DimType typeIn, std::string speciesIn);
		Dim(std::string nameIn, double minIn, double maxIn, int nIn, DimType typeIn, std::string species1In, std::string species2In);
		Dim(const Dim& d);
		Dim & operator=(const Dim& d);
		~Dim();

		// Copy function
		void _copy(const Dim& d);

		// Dimensions to use as the basis function
		bool bf_use_all_dim; // By default, use all available
		// Otherwise, use these:
		int bf_n_dim;
		std::vector<int> bf_dim_idxs;
		std::vector<std::string> bf_dim_names;
		// Specify dimensions
		void set_basis_func_dims(std::string dim_name);
		void set_basis_func_dims(std::vector<std::string> dim_name);

		// Update moments from species
		enum MomentType {AWAKE, ASLEEP};
		void append_moments_from_latt(MomentType moment_type);

		// Comparators
		bool operator==(const Dim& d);
		bool operator!=(const Dim& d);
	};

	/****************************************
	n-Dimensional grid
	****************************************/

	class GridND
	{
	private:
		
		// Internal copy func
		void _copy(const GridND& g);

	protected:

		// Number of dimensions
		int _n_dim;

		// Vals
		double *_vals;
		int _val_len;

		// Dimensions
		std::vector<Dim> _dims;

		// Dimension length squares
		std::vector<int> _dim_pwrs;

	public:

		// Constructor
		GridND(std::string name, int n_dim, std::vector<Dim> dims);
		GridND(std::string name, Dim dim);
		GridND(const GridND& bf);
		GridND& operator=(const GridND& bf);
		~GridND();

		// Name - useful for writing!
		std::string name;

		// Get/set an element by index
		double get(int *idxs) const;
		double get(int i) const;
		void set(int *idxs, double val);
		void set(int i, double val);

		// Multiply all values by some constant
		void multiply_all(double val);

		// Increment by a whole grid
		void increment(const GridND& g);

		// Get indexes by element
		int* get_idxs(int i);

		// Write to a file
		void write_to_file(std::string dir, int idx);

		// Check dimensions against another grid
		bool check_dims(const GridND& g);

		// Zero
		void zero();
	};

	/****************************************
	BasisFuncND
	****************************************/

	class BasisFuncND : public GridND {

	private:

		// Get bounding n-dim cube of 4 pts
		void _get_bounding(double *x);
		void _fill_p(int dim);
		// Parameters - only alloc/dealloc only once
		int* _idxs_bounding;
		int* _idxs_p_cube;
		int* _idxs_ext_1;
		int* _idxs_ext_2;
		double *_fracs;
		double *_p_cube;

		// Internal copy/clean up function
		void _copy(const BasisFuncND& bf);
		void _clean_up();

		// Interpolation
		double _cubic_interp(double p[4], double f);
		double _deriv(double p[4], double f);
		double _n_cubic_interp(int dim, double* p, double fracs[], bool derivs[]);

	public:

		// Constructor
		BasisFuncND(std::string name, int n_dim, std::vector<Dim> dims);
		BasisFuncND(const BasisFuncND& bf);
		BasisFuncND& operator=(const BasisFuncND& bf);
		~BasisFuncND();

		// Get values, if they are in the lattice
		double get_val(double *x);
		double get_deriv(double *x, bool *derivs);

		// Update with delta f
		void update(const GridND &g);

		// Test fill in various dimensions
		double test_func_2d(double x, double y);
		void test_fill_2d();
		double test_func_3d(double x, double y, double z);
		void test_fill_3d();

		// From parent
		
		// Get/set an element by index
		double get(int *idxs) const;
		double get(int i) const;
		void set(int *idxs, double val);

		// Write to a file
		void write_to_file(std::string dir, int idx);

	};
};

