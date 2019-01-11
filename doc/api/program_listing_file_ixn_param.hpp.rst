
.. _program_listing_file_ixn_param.hpp:

Program Listing for File ixn_param.hpp
======================================

- Return to documentation for :ref:`file_ixn_param.hpp`

.. code-block:: cpp

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
       Grid
       ****************************************/
   
       class Grid {
   
       private:
   
           // Name
           std::string _name;
   
           // Size
           int _n;
           // Min/max/increment
           double _min;
           double _max;
           double _delta;
           // Grid of vals
           double *_grid;
   
           // Copy/clean up
           void _copy(const Grid& other);
           void _clean_up();
   
       public:
   
           // Constructors
           Grid(std::string name, double min, double max, int n);
           Grid(const Grid& other);
           Grid & operator=(const Grid& other);
           ~Grid();
   
           // Getters/setters
           std::string name() const;
           double delta() const;
           int n() const;
           double get_by_idx(int i) const;
   
           // Get indexes surrounding a point
           // ie point is between i and i+1 where i is returned
           int surrounding_idxs(double x) const; 
   
           // Get fraction of a point between successive points
           double frac_between(double x) const;
           // Second optional specification: the return of the surrounding idxs
           double frac_between(double x, int i) const;
   
           // Check if a given point is in the grid
           bool in_grid(double x) const;
   
           // Print grid range
           void print_grid_range() const;
   
           // Write the grid into an ofstream
           void write_grid(std::string fname) const;
   
           // Test: values if a sin function were defined on the grid
           std::vector<double> test_sin() const;
           std::vector<double> test_cos() const;
       };
   
       /****************************************
       Forward declare
       ****************************************/
   
       class BasisFunc;
       class VarTerm;
   
       /****************************************
       Interaction parameter
       ****************************************/
   
       // Enumeration of type of dimension
       enum IxnParamType { Hp, Jp };
   
       class IxnParam : public Grid {
   
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
           void _copy(const IxnParam& other);
           void _clean_up();
   
       public:
   
           // Constructor
           IxnParam(std::string name, IxnParamType type, Species *sp, double min, double max, int n, double val0, int n_t);
           IxnParam(std::string name, IxnParamType type, Species *sp1, Species *sp2, double min, double max, int n, double val0, int n_t);
           IxnParam(const IxnParam& other);
           IxnParam & operator=(const IxnParam& other);
           ~IxnParam();
   
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
   
       /****************************************
       Array
       ****************************************/
   
       class Array {
       private:
           
           // Internal copy func/clean up
           void _copy(const Array& other);
           void _clean_up();
   
       protected:
   
           // Number of dimensions
           int _n_params;
   
           // Vals
           double *_vals;
           int _val_len;
   
           // Dimensions
           std::vector<IxnParam*> _ixn_params;
   
           // Dimension length squares
           std::vector<int> _dim_pwrs;
   
       public:
   
           // Constructor
           Array(std::vector<IxnParam*> ixn_params);
           Array(IxnParam* ixn_param);
           Array(const Array& other);
           Array& operator=(const Array& other);
           Array(Array&& other);
           Array& operator=(Array&& other);
           ~Array();
   
           // Get/set an element by index
           double get_by_idxs(int *idxs) const;
           double get_by_idx(int i) const;
           void set_by_idxs(int *idxs, double val);
           void set_by_idx(int i, double val);
   
           // Get indexes by element
           void get_idxs(int i, int* idxs) const;
   
           // Write/Read to a file
           void write_grid(std::string fname) const;
           void write_vals(std::string dir, std::string name, int idx) const;
           void read_vals(std::string fname);
   
           // Check dimensions against another array
           bool check_dims(const Array& other) const;
   
           // Zero
           void zero();
       };
   
       /****************************************
       Variational Term Trajectory
       ****************************************/
   
       class VarTerm {
       private:
   
           // Name
           std::string _name;
   
           // Length of trajectory
           int _n_t;
   
           // Numerator and denominator
           IxnParam *_num;
           BasisFunc *_denom;
   
           // Whether or not this term deserves a delta source
           bool _delta_source;
   
           // Basis func corresponding to num
           BasisFunc *_num_bf;
   
           // Derivs of the numerators basis function necessary for updating
           int _n_ixn_params_in_num_bf;
           double *_num_bf_derivs;
   
           // Values over time
           std::vector<Array> _vals;
   
           // Val length of each array
           int _val_len;
   
           // Pointers to other ixn params and variational terms needed to update this term
           // The interaction parameters that are arguments to _num_bf, so we can take derivatives
           std::vector<VarTerm*> _update_var_terms;
   
           // Clean up/copy
           void _copy(const VarTerm& other);
           void _clean_up();
   
       public:
   
           // Constructor
           VarTerm(std::string name, IxnParam *num, BasisFunc *denom, std::vector<IxnParam*> denom_ixn_params, BasisFunc *num_bf, int n_ixn_params_in_num_bf, int n_t);
           VarTerm(const VarTerm& other);
           VarTerm & operator=(const VarTerm& other);
           ~VarTerm();
   
           // Set the pointers needed to update this term
           void add_update_ptr(VarTerm* var_term);
   
           // Validate setup
           void validate_setup() const;
   
           // Calculate next timestep
           void calculate_at_time(int it_next, double dt);
   
           // Getters/setters
           double get_at_time_by_idx(int it, int i);
           std::string name();
   
           // Write
           void write_vals(std::string dir, int idx) const;
       };
   
       /****************************************
       BasisFunc
       ****************************************/
   
       class BasisFunc : public Array {
   
       private:
   
           // Name
           std::string _name;
   
           // The variational terms and ixn_params needed to update
           std::vector<std::pair<IxnParam*,VarTerm*>> _update_ptrs;
   
           // Get bounding n-dim cube of 4 pts
           void _get_bounding(int it, bool safe=false);
           void _fill_p(int dim);
           // Parameters - only alloc/dealloc only once
           int* _idxs_bounding;
           int* _idxs_p_cube;
           int* _idxs_ext_1;
           int* _idxs_ext_2;
           double *_fracs;
           double *_p_cube;
   
           // Update, if needed
           double *_update_gathered;
   
           // Derivatives
           bool *_derivs;
   
           // Internal copy/clean up function
           void _copy(const BasisFunc& other);
           void _clean_up();
   
           // Interpolation
           double _cubic_interp(double p[4], double f);
           double _deriv(double p[4], double f);
           double _n_cubic_interp(int dim, double* p, double fracs[], bool derivs[]);
   
       public:
   
           // Constructor
           BasisFunc(std::string name, std::vector<IxnParam*> ixn_params);
           BasisFunc(const BasisFunc& other);
           BasisFunc& operator=(const BasisFunc& other);
           ~BasisFunc();
   
           // Add pointers needed to update
           void add_update_ptrs(IxnParam* ixn_param, VarTerm* var_term);
   
           // Validate setup
           void validate_setup() const;
   
           // Get values, if they are in the lattice
           double get_at_time(int it);
           double get_deriv_at_time(int it, int i_dim);
   
           // Name
           std::string name() const;
   
           // Calculate the new basis function
           void update(int n_t, double dt, double dopt);
           void update_gather(int n_t, double dt, double dopt);
           void update_committ_gathered();
   
           // Test fill in various dimensions
           void test_fill_2d();
           void test_fill_3d();
   
           // From parent
           
           // Get/set an element by index
           double get_by_idxs(int *idxs) const;
           double get_by_idx(int i) const;
           void set_by_idxs(int *idxs, double val);
   
           // Write/Read grid/vals
           void write_grid(std::string fname) const;
           void write_vals(std::string dir, int idx) const;
           void read_vals(std::string fname);
   
           // Get the delta source
           double get_delta_source(int it, int i);
       };
   };
   
