#ifndef VECTOR_h
#define VECTOR_h
#include <vector>
#endif

#ifndef STRING_h
#define STRING_h
#include <string>
#endif

#ifndef HIDDEN_UNIT_h
#define HIDDEN_UNIT_h
#include "hidden_unit.hpp"
#endif

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Interaction parameter
	****************************************/

	// Enumeration of type of dimension
	enum IxnParamType { Hp, Jp, Wp, Kp, Bp };

	class IxnParam {

	private:

		// Name
		std::string _name;

		// Type
		IxnParamType _type;

		// Species
		// If Hp or Wp: only sp1
		Species *_sp1;
		// If Jp or Kp: also sp2
		Species *_sp2;
		// If Kp: also sp3
		Species *_sp3;
		// All species
		bool _all_species;
		std::vector<Species*> _sp_all;

		// If Wp:
		// Connects sites (visible) to hidden units
		std::vector<std::pair<Site*,HiddenUnit*> > _conns;

		// If Bp:
		// List of hidden units
		std::vector<HiddenUnit*> _hidden_units;

		// Number of time points in these trajs
		int _n_t;

		// Value
		double _val;

		// Solution traj
		bool _track_soln_traj;
		std::vector<double> _soln_traj;

		// Initial guess
		double _val_guess;

		// Awake and asleep moments
		double _asleep;
		double _awake;

		// Copy, clean up
		void _copy(const IxnParam& other);
		void _reset();
		void _clean_up();

	public:

		// Constructor
		IxnParam(std::string name, IxnParamType type, std::vector<Species*> sp_all, double val_guess);
		IxnParam(std::string name, IxnParamType type, Species *sp, double val_guess);
		IxnParam(std::string name, IxnParamType type, Species *sp1, Species *sp2, double val_guess);
		IxnParam(std::string name, IxnParamType type, Species *sp1, Species *sp2, Species *sp3, double val_guess);
		IxnParam(const IxnParam& other);
		IxnParam(IxnParam&& other);
		IxnParam & operator=(const IxnParam& other);
		IxnParam & operator=(IxnParam&& other);
		~IxnParam();

		// Check if this ixn param is...
		bool is_type_with_species(IxnParamType type, std::string s) const;
		bool is_type_with_species(IxnParamType type, std::string s1, std::string s2) const;
		bool is_type_with_species(IxnParamType type, std::string s1, std::string s2, std::string s3) const;
		bool is_type_for_any_species(IxnParamType type) const;

		// If Wp
		// Add a visible->hidden unit connection
		void add_visible_hidden_connection(Site *sptr, HiddenUnit *hup);

		// If Bp
		// Add a hidden unit to monitor
		void add_hidden_unit(HiddenUnit *hup);

		// Update based on diff
		void update(double dopt, bool l2_reg=false, double lambda=0.);

		// Getters/setters
		std::string name() const;
		double get() const;
		void set_val(double val);
		void set_guess(double guess);
		void reset();
		std::vector<Species*> get_species() const;

		// Moments from lattice
		enum MomentType {AWAKE, ASLEEP};
		double get_moment(MomentType moment_type) const;
		void moments_reset(MomentType moment_type);
		void moments_retrieve(MomentType moment_type);
		void moments_retrieve(MomentType moment_type, int batch_size);
		double moments_diff() const;

		// Store solution traj
		void set_track_soln_traj(bool flag);
		void reset_soln_traj();
		double get_ave() const;
		double get_ave(int last_n_steps) const;
	};
};

