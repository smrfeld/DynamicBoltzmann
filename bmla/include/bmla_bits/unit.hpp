#include <unordered_map>

#ifndef FWDS_SPECIES_H
#define FWDS_SPECIES_H
#include "fwds/fwds_species.hpp"
#endif

/************************************
* Namespace for bmla
************************************/

namespace bmla {

	// Forward
	class BiasDict;
	class ConnVV;
	class ConnVVV;
	class ConnVH;

	/****************************************
	Class to hold a lattice site
	--- CAUTION : abstract base! ---
	***************************************/

	class Unit {

	private:

		class Impl;
		std::unique_ptr<Impl> _impl;		

	public:

		/********************
		Constructor
		********************/

		Unit(int x);
		Unit(int x, int y);
		Unit(int x, int y, int z);
		Unit(int x, std::vector<Sptr> species_possible);
		Unit(int x, int y, std::vector<Sptr> species_possible);
		Unit(int x, int y, int z, std::vector<Sptr> species_possible);
		Unit(const Unit& other);
		Unit(Unit&& other);
		Unit& operator=(const Unit& other);
		Unit& operator=(Unit&& other);
		~Unit();

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
		bool less_than(const Unit &other) const;

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
		--- CAUTION: pure virtual ---
		********************/

		virtual double get_activation_for_species(Sptr &species) const = 0;

		/********************
		Sample
		--- CAUTION: pure virtual ---
		********************/

		virtual void sample(bool binary=true) = 0;

		// Sample given activation for every species
		void sample_given_activations(const std::vector<double>& activations, bool binary=true);
	};
	// Comparator
	bool operator <(const Unit& a, const Unit& b);


































	/****************************************
	Visible
	***************************************/

	class UnitVisible : public Unit {

	private:

		class Impl;
		std::unique_ptr<Impl> _impl;		

	public:

		/********************
		Constructor
		********************/

		// using Unit::Unit;
		UnitVisible(int x);
		UnitVisible(int x, int y);
		UnitVisible(int x, int y, int z);
		UnitVisible(int x, std::vector<Sptr> species_possible);
		UnitVisible(int x, int y, std::vector<Sptr> species_possible);
		UnitVisible(int x, int y, int z, std::vector<Sptr> species_possible);
		UnitVisible(const UnitVisible& other);
		UnitVisible(UnitVisible&& other);
		UnitVisible& operator=(const UnitVisible& other);
		UnitVisible& operator=(UnitVisible&& other);
		~UnitVisible();

		/********************
		Add connection
		********************/

		void add_conn(ConnVV *conn, int idx_of_me); // idx = 0,1
		const std::vector<std::pair<ConnVV*,int>>& get_conns_vv() const;

		void add_conn(ConnVVV *conn, int idx_of_me); // idx = 0,1,2
		const std::vector<std::pair<ConnVVV*,int>>& get_conns_vvv() const;
		
		void add_conn(ConnVH *conn);
		const std::vector<ConnVH*>& get_conns_vh() const;

		/********************
		Get activation
		********************/

		double get_activation_for_species(Sptr &species) const;

		/********************
		Sample
		********************/

		void sample(bool binary=true);

	};

































	/****************************************
	Hidden
	***************************************/

	class UnitHidden : public Unit {

	private:

		class Impl;
		std::unique_ptr<Impl> _impl;		

	public:

		/********************
		Constructor
		********************/

		UnitHidden(int layer, int x);
		UnitHidden(int layer, int x, int y);
		UnitHidden(int layer, int x, int y, int z);
		UnitHidden(int layer, int x, std::vector<Sptr> species_possible);
		UnitHidden(int layer, int x, int y, std::vector<Sptr> species_possible);
		UnitHidden(int layer, int x, int y, int z, std::vector<Sptr> species_possible);
		UnitHidden(const UnitHidden& other);
		UnitHidden(UnitHidden&& other);
		UnitHidden& operator=(const UnitHidden& other);
		UnitHidden& operator=(UnitHidden&& other);
		~UnitHidden();

		/********************
		Location
		********************/

		int layer() const;

		/********************
		Finish setup in lattice
		********************/

		void add_conn(ConnVH *conn);
		const std::vector<ConnVH*>& get_conns_vh() const;

		/********************
		Get activation
		********************/

		double get_activation_for_species(Sptr &species) const;

		/********************
		Sample
		********************/

		void sample(bool binary=true);
	};


};