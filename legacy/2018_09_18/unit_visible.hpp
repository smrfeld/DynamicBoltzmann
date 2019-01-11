#ifndef UNORDERED_MAP_H
#define UNORDERED_MAP_H
#include <unordered_map>
#endif

#ifndef FWDS_SPECIES_H
#define FWDS_SPECIES_H
#include "fwds/fwds_species.hpp"
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	// Forward
	class BiasDict;
	class ConnVV;
	class ConnVVV;
	class ConnVH;

	/****************************************
	Class to hold a lattice site
	***************************************/

	class UnitVisible {

	private:

		class Impl;
		std::unique_ptr<Impl> _impl;		

	public:

		/********************
		Constructor
		********************/

		UnitVisible(int x);
		UnitVisible(int x, std::vector<Sptr> species_possible);
		UnitVisible(int x, int y);
		UnitVisible(int x, int y, std::vector<Sptr> species_possible);
		UnitVisible(int x, int y, int z);
		UnitVisible(int x, int y, int z, std::vector<Sptr> species_possible);
		UnitVisible(const UnitVisible& other);
		UnitVisible(UnitVisible&& other);
		UnitVisible& operator=(const UnitVisible& other);
		UnitVisible& operator=(UnitVisible&& other);
		~UnitVisible();

		/********************
		Check setup
		********************/

		void check_setup() const;

		/********************
		Location
		********************/

		int dim() const;
		int x() const;
		int y() const;
		int z() const;
		bool less_than(const UnitVisible &other) const;

		/********************
		Finish setup in lattice
		********************/

		// Add possible species
		void add_possible_species(Sptr species);
		void set_possible_species(std::vector<Sptr> species);
		const std::vector<Sptr>& get_possible_species() const;

		// Bias dict
		void set_bias_dict(std::shared_ptr<BiasDict> bias_dict);
		const std::shared_ptr<BiasDict>& get_bias_dict() const;

		// Add connection
		void add_conn(ConnVV *conn, int idx_of_me); // idx = 0,1
		const std::vector<std::pair<ConnVV*,int>>& get_conns_vv() const;

		void add_conn(ConnVVV *conn, int idx_of_me); // idx = 0,1,2
		const std::vector<std::pair<ConnVVV*,int>>& get_conns_vvv() const;
		
		void add_conn(ConnVH *conn);
		const std::vector<ConnVH*>& get_conns_vh() const;

		/********************
		Get/set probability of different species to live on this site
		********************/

		// Check mode
		bool check_is_b_mode() const;

		// Flip between the two modes
		void set_b_mode(bool flag);

		// Binary
		Sptr get_b_mode_species() const; // nullptr for empty
		void set_b_mode_species(Sptr sp);
		void set_b_mode_species(std::string sp);
		void set_b_mode_empty();
		bool check_b_mode_is_empty() const;

		// Probabilistic
		double get_p_mode_prob(Sptr sp) const; // nullptr for empty
		const std::unordered_map<Sptr, double>& get_p_mode_probs() const;
		void set_p_mode_prob(Sptr sp, double prob);
		void set_p_mode_prob(std::string sp, double prob);
		void set_p_mode_empty();
		bool check_p_mode_is_empty() const;

		// Convert between the two
		void convert_b_to_p_mode();
		void convert_p_to_b_mode();

		/********************
		Get moment
		********************/

		double get_moment(std::string ixn_param_name, bool binary=true) const;

		/********************
		Get activation
		********************/

		double get_act_for_species_at_timepoint(Sptr &species, int timepoint) const;

		/********************
		Sample
		********************/

		void sample_at_timepoint(int timepoint, bool binary=true);
	};
	// Comparator
	bool operator <(const UnitVisible& a, const UnitVisible& b);

};