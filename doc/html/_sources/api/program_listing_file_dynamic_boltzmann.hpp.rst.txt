
.. _program_listing_file_dynamic_boltzmann.hpp:

Program Listing for File dynamic_boltzmann.hpp
==============================================

- Return to documentation for :ref:`file_dynamic_boltzmann.hpp`

.. code-block:: cpp

   #ifndef IXN_PARAM_H
   #define IXN_PARAM_H
   #include "ixn_param.hpp"
   #endif
   
   #ifndef LATTICE_H
   #define LATTICE_H
   #include "lattice.hpp"
   #endif
   
   #ifndef LIST_H
   #define LIST_H
   #include <list>
   #endif
   
   #define DIAG_SETUP 0
   #define DIAG_SOLVE 0
   
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
   
           // Min/max/npts
           double min,max;
           int n;
   
           // Initial value
           double init;
   
           // Basis function dimension names
           std::vector<std::string> basis_func_dims;
   
           // Constructor
           Dim(std::string name, DimType type, std::string species, std::vector<std::string> basis_func_dims, double min, double max, int n, double init);
           Dim(std::string name, DimType type, std::string species1, std::string species2, std::vector<std::string> basis_func_dims, double min, double max, int n, double init);
       };
   
       /****************************************
       OptProblem
       ****************************************/
   
       class OptProblem {
   
       private:
   
           /********************
           Parameters
           ********************/
   
           // The directory to write to
           std::string _dir_io;
   
           // Number of dimensions
           int _n_param;
   
           // List of interaction parameters
           std::list<IxnParam> _ixn_params;
   
           // List of basis funcs
           std::list<BasisFunc> _bfs;
   
           // List of variational terms
           std::list<VarTerm> _var_terms;
   
           // Time dimension
           Grid _time;
   
           // Species present
           std::list<Species> _species;
   
           // Filenames to choose from
           std::vector<std::string> _fnames;
   
           // Number of steps in this nu solution
           int _n_t_soln;
   
           // The current time in the optimization step
           int _t_opt;
   
           // Batch size
           int _n_batch;
   
           // Number of annealing steps
           int _n_annealing;
   
           // Lattice size
           int _box_length;
   
           // Lattice to hold the current sample of the batch
           Lattice _latt;
   
           // Update step for optimization
           double _dopt;
   
           // Number opt steps
           int _n_opt;
   
           // Start index for the files to read
           int _fname_start_idx;
   
           // Flags
           bool _write_bf_only_last;
   
           /********************
           Search functions
           ********************/
   
           Species* _find_species(std::string name);
           IxnParam* _find_ixn_param(std::string name);
           BasisFunc* _find_basis_func(std::string name);
           VarTerm* _find_var_term(std::string name);
   
           /********************
           Helpers
           ********************/
   
           void _clean_up();
           void _copy(const OptProblem& other);
   
       public:
   
           /********************
           Constructor
           ********************/
   
           OptProblem(std::vector<Dim> dims, std::vector<std::string> species, double t_max, int n_t, int batch_size, int n_annealing, int box_length, double dopt, int n_opt);
   
           OptProblem(const OptProblem& other);
   
           OptProblem(OptProblem&& other);
   
           OptProblem & operator=(const OptProblem& other);
   
           OptProblem& operator=(OptProblem&& other)
   
           
           ~OptProblem();
   
           /********************
           Set properties  
           ********************/
   
           void set_dir_io(std::string dir);
   
           void set_fname_start_idx(int idx);
           
           void add_fname(std::string f);
   
           /********************
           Validate setup
           ********************/
   
           void validate_setup() const;
   
           /********************
           Solve interaction parameter traj
           ********************/
   
           void solve_ixn_param_traj();
   
           /********************
           Solve for variational trajectory
           ********************/
   
           void solve_var_traj();
   
           /********************
           Solve
           ********************/
   
           void solve(bool verbose=false);
           void solve_varying_ic(bool verbose=false);
   
           /********************
           Read some initial conditions
           ********************/
   
           void read_init_cond(std::string dir);
   
           /********************
           Write
           ********************/
   
           void write_bf_grids() const;
           void write_t_grid() const;
   
           void write_ixn_params(std::string dir, int idx) const;
           void write_ixn_params(std::string dir, int idx1, int idx2) const;
           void write_bfs(std::string dir, int idx) const;
           void write_var_terms(std::string dir, int idx) const;
           void write_moments(std::string dir, int idx) const;
           void write_moments(std::string dir, int idx1, int idx2) const;
   
           void set_flag_write_bf_only_final();
   
           /********************
           Read
           ********************/
   
           void read_bf(std::string bf_name, std::string fname);
       };
   
   
   };
