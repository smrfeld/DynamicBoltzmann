#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

#ifndef STRING_H
#define STRING_H
#include <string>
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	Interaction parameter
	****************************************/

	// Enumeration of type of dimension
	enum IxnParamType { Hp, Jp, Kp, Wp, Bp };

	// Forward
	class Species;
	struct Species2;
	struct Species3;
	struct SpeciesVH;
	class HiddenSpecies;
	class Site;
	class HiddenUnit;
	class BasisFunc;
	class Domain1D;

	class IxnParamTraj {

	private:

		// Name
		std::string _name;

		// Type
		IxnParamType _type;

		// Domain
		double _domain_min;
		double _domain_max;
		double _domain_delta;
		int _domain_no_pts;

		// Species
		
		// For a visible bias
		std::vector<Species*> _sp_bias_visible;
		// For a hidden bias
		std::vector<HiddenSpecies*> _sp_bias_hidden;
		// For a visible J
		std::vector<Species2> _sp_doublet;
		// For a visible K
		std::vector<Species3> _sp_triplet;
		// For a visible to hidden weight W
		std::vector<SpeciesVH> _sp_conn;

		// If Wp:
		// Connects sites (visible) to hidden units
		std::vector<std::pair<Site*,HiddenUnit*> > _conns;

		// If Bp:
		// List of hidden units
		std::vector<HiddenUnit*> _hidden_units;

		// Number of time points in these trajs
		int _n_t;

		// Values over time
		double *_vals;

		// Initial value
		double _val0;

		// Awake and asleep moments over time
		double *_asleep;
		double *_awake;

		// Fixed value for the awake moment
		bool _is_awake_fixed;
		double *_awake_fixed;

		// Pointer to the basis function
		BasisFunc *_bf;

		// Pointer to the optimization time
		int *_t_opt_ptr;

		// Copy, clean up
		void _clean_up();
		void _reset();
		void _copy(const IxnParamTraj& other);

	public:

		/********************
		Constructor
		********************/

		IxnParamTraj(std::string name, IxnParamType type, double domain_min, double domain_max, int domain_no_pts, double val0, int n_t, int *t_opt_ptr);
		IxnParamTraj(const IxnParamTraj& other);
		IxnParamTraj(IxnParamTraj&& other);
		IxnParamTraj& operator=(const IxnParamTraj& other);
		IxnParamTraj& operator=(IxnParamTraj&& other);
		~IxnParamTraj();

		/********************
		Name
		********************/

		std::string get_name() const;

		/********************
		Domain
		********************/

		// Getters
		int get_domain_no_pts() const;
		double get_domain_min() const;
		double get_domain_max() const;
		double get_domain_delta() const;

		// Get pt in domain
		double get_domain_pt_by_idx(int i) const;

		// Check if point is in domain
		bool check_if_pt_is_inside_domain(double x) const;

		// Get indexes surrounding a point
		// ie point is between i and i+1 where i is returned
		int get_domain_idxs_surrounding_pt(double x) const; 

		// Get fraction of a point between successive points
		double get_domain_frac_between(double x) const;
		// Second optional specification: the return of the surrounding idxs
		double get_domain_frac_between(double x, int i) const;

		// Print domain range
		void print_domain_range() const;

		/********************
		Set time
		********************/

		void set_n_t(int n_t);

		/********************
		Add a species
		********************/

		void add_species(Species *sp);
		void add_species(HiddenSpecies *hsp);
		void add_species(Species *sp1, Species *sp2);
		void add_species(Species *sp1, Species *sp2, Species *sp3);
		void add_species(Species *sp, HiddenSpecies *hsp);

		/********************
		If W: Add a visible->hidden unit connection
		********************/

		void add_visible_hidden_connection(Site *sptr, HiddenUnit *hup);

		/********************
		If B: Add a hidden unit to monitor
		********************/

		void add_hidden_unit(HiddenUnit *hup);

		/********************
		Set/check basis func pointer
		********************/

		void set_basis_func_ptr(BasisFunc* bf);
		bool is_bf(BasisFunc *bf);

		/********************
		Set IC
		********************/

		void set_init_cond(double val);

		/********************
		Set fixed awake
		********************/

		void set_fixed_awake_moment(std::vector<double> vals);

		/********************
		Validate setup
		********************/

		void validate_setup() const;

		/********************
		Getters
		********************/

		// At the current time in the optimization
		double get() const;
		// At a specific time
		double get_at_time(int it) const;

		/********************
		Get species involved
		********************/

		std::vector<Species*> get_species_bias() const;
		std::vector<HiddenSpecies*> get_species_hidden_bias() const;
		std::vector<Species2> get_species_doublet() const;
		std::vector<Species3> get_species_triplet() const;
		std::vector<SpeciesVH> get_species_conn() const;

		/********************
		Calculate the next step
		********************/

		// Returns false if new point is out of grid
		bool calculate_at_time(int it_next, double dt);

		/********************
		Moments from lattice
		********************/

		enum MomentType {AWAKE, ASLEEP};
		void moments_reset();
		void moments_retrieve_at_time(MomentType moment_type, int it);
		void moments_retrieve_at_time(MomentType moment_type, int it, int batch_size);
		double moments_diff_at_time(int it);
		double moments_get_at_time(MomentType moment_type, int it) const;

		/********************
		Write into an ofstream
		********************/

		void write_vals(std::string dir, int idx_opt_step, std::vector<int> idxs, int n_t_traj) const;
		void write_moments(std::string dir, int idx_opt_step, std::vector<int> idxs, int n_t_traj) const;
	};
};

