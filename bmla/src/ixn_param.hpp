#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

#ifndef STRING_H
#define STRING_H
#include <string>
#endif

/************************************
* Namespace for bmla
************************************/

namespace bmla {

	// Forwards
	class Species;
	class HiddenSpecies;
	struct Species2;
	struct Species3;
	struct Species4;
	struct SpeciesVH;
	class Site;
	class HiddenUnit;

	/****************************************
	Interaction parameter
	****************************************/

	// Enumeration of type of dimension
	enum IxnParamType { Hp, Jp, Wp, Kp, Bp, Qp };

	class IxnParam {

	private:

		// Name
		std::string _name;

		// Type
		IxnParamType _type;

		// Species
		
		// For a visible bias
		std::vector<Species*> _sp_bias_visible;
		// For a hidden bias
		std::vector<HiddenSpecies*> _sp_bias_hidden;
		// For a visible J
		std::vector<Species2> _sp_doublet;
		// For a visible K
		std::vector<Species3> _sp_triplet;
		// For a visible Q
		std::vector<Species4> _sp_quartic;
		// For a visible to hidden weight W
		std::vector<SpeciesVH> _sp_conn;

		// If Wp:
		// Connects sites (visible) to hidden units
		std::vector<std::pair<Site*,HiddenUnit*> > _conns;

		// If Bp:
		// List of hidden units
		std::vector<HiddenUnit*> _hidden_units;

		// Value
		double _val;

		// Nesterov: previous value
		double _val_nesterov;

		// Solution traj
		bool _track_soln_traj;
		std::vector<double> _soln_traj;

		// Initial guess
		double _val_guess;

		// Awake and asleep moments
		double _asleep;
		double _awake;

		// Fixed value for the awake moment
		bool _is_awake_fixed;
		double _awake_fixed;

		// Copy, clean up
		void _copy(const IxnParam& other);
		void _reset();
		void _clean_up();

	public:

		/********************
		Constructor
		********************/

		IxnParam(std::string name, IxnParamType type, double val_guess);
		IxnParam(const IxnParam& other);
		IxnParam(IxnParam&& other);
		IxnParam & operator=(const IxnParam& other);
		IxnParam & operator=(IxnParam&& other);
		~IxnParam();

		/********************
		Nesterov
		********************/

		// Move to the nesterov intermediate point
		void nesterov_move_to_intermediate_pt(int opt_step);

		// Set prev nesterov
		void nesterov_set_prev_equal_curr();

		/********************
		Add a species
		********************/

		void add_species(Species *sp);
		void add_species(HiddenSpecies *hsp);
		void add_species(Species *sp1, Species *sp2);
		void add_species(Species *sp1, Species *sp2, Species *sp3);
		void add_species(Species *sp1, Species *sp2, Species *sp3, Species *sp4);
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
		Update based on diff
		********************/

		void update(double dopt, bool l2_reg=false, double lambda=0.);

		/********************
		Setters/getters
		********************/

		std::string name() const;
		double get() const;
		void set_val(double val);
		void set_guess(double guess);
		void set_fixed_awake_moment(double val);
		void reset_to_guess();

		/********************
		Get species involved
		********************/

		std::vector<Species*> get_species_bias() const;
		std::vector<HiddenSpecies*> get_species_hidden_bias() const;
		std::vector<Species2> get_species_doublet() const;
		std::vector<Species3> get_species_triplet() const;
		std::vector<Species4> get_species_quartic() const;
		std::vector<SpeciesVH> get_species_conn() const;

		/********************
		Moments from lattice
		********************/

		enum MomentType {AWAKE, ASLEEP};
		double get_moment(MomentType moment_type) const;
		void moments_reset(MomentType moment_type);
		void moments_retrieve(MomentType moment_type);
		void moments_retrieve(MomentType moment_type, int batch_size);
		double moments_diff() const;

		/********************
		Store solution traj
		********************/

		void set_track_soln_traj(bool flag);
		void reset_soln_traj();
		double get_ave() const;
		double get_ave(int last_n_steps) const;
	};
};

