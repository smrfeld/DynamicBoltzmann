#include <unordered_map>
#include <map>
#include <vector>

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
	class ConnHH;

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

		virtual void print() const;
		virtual std::string print_str() const;

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

		double get_occ(Sptr sp) const; // nullptr for empty
		const std::unordered_map<Sptr, double>& get_nonzero_occs() const;
		void set_occ(Sptr sp, double prob);
		void set_occ(std::string sp, double prob);
		void set_occ_random();

		bool get_is_empty() const;
		void set_empty();

		/********************
		Get moment
		********************/

		double get_moment(std::string ixn_param_name) const;

		/********************
		Sample
		********************/

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
		Verbose
		********************/

		void print() const;
		std::string print_str() const;

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

		double get_activation_for_species_at_timepoint(Sptr &species, int timepoint) const;

		/********************
		Sample
		********************/

		void sample_at_timepoint(int timepoint, bool binary);

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

		UnitHidden(int layer, int idx);
		UnitHidden(int layer, int idx, std::vector<Sptr> species_possible);
		UnitHidden(const UnitHidden& other);
		UnitHidden(UnitHidden&& other);
		UnitHidden& operator=(const UnitHidden& other);
		UnitHidden& operator=(UnitHidden&& other);
		~UnitHidden();

		/********************
		Verbose
		********************/

		void print() const;
		std::string print_str() const;

		/********************
		Location
		********************/

		int layer() const;
		int idx() const;

		/********************
		Finish setup in lattice
		********************/

		// Visible-hidden conns
		void add_conn(ConnVH *conn);
		const std::vector<ConnVH*>& get_conns_vh() const;

		// Hidden-hidden conns
		void add_conn(ConnHH *conn, int idx_of_me, int from_layer);
		const std::map<int,std::vector<std::pair<ConnHH*,int>>>& get_conns_hh() const;

		/********************
		Get activation
		********************/

		double get_activation_for_species_at_timepoint(Sptr &species, int timepoint, int given_layer) const;

		/********************
		Sample
		********************/

		void sample_at_timepoint(int timepoint, bool binary);
		void sample_at_timepoint(int timepoint, bool binary, int given_layer);
	};
};
