#include <unordered_map>
#include <map>
#include <vector>

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
	class ConnHH;

	/****************************************
	Class to hold a lattice site
	***************************************/

	struct Act {
		// Species
		// nullptr for empty
		Sptr sp;
		// activation
		double act;
		// prob = exp(activation)
		double prob;
		// prop = cumulative sum
		double prop;
	};

	class Unit {

	private:

		// Dimensionality and location
		int _dim;
		int _x,_y,_z;

		// Probabilistic mode
		std::unordered_map<Sptr, double> _nonzero_occs, _nonzero_occs_tbc;

		// Data structures for sampling
		double _sampling_rand;
		int _sampling_i_chosen;

		// Sample a vector of propensities (cumulative probabilities)
		void _sample_prop_vec();
		void _prepare_sample(bool binary);

		// Constructor helpers
		void _clean_up();
		void _move(Unit& other);
		void _copy(const Unit& other);

	protected:

		// Activations
		std::vector<Act> _activations;

		// Bias dict
		std::shared_ptr<BiasDict> _bias_dict;

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

		void print() const;
		std::string print_str() const;

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

		// Bias dict
		void set_bias_dict(std::shared_ptr<BiasDict> bias_dict);
		const std::shared_ptr<BiasDict>& get_bias_dict() const;

		/********************
		Get/set probability of different species to live on this site
		********************/

		double get_occ(Sptr sp) const; // nullptr for empty
		const std::unordered_map<Sptr, double>& get_nonzero_occs() const;
		void set_occ(Sptr sp, double prob);
		void set_occ_random();

		void binarize();

		bool get_is_empty() const;
		void set_empty();

		/********************
		Get moment
		********************/

		double get_moment(std::string ixn_param_name) const;

		/********************
		Sample
		********************/

		virtual void form_propensity_vector();
		virtual void form_propensity_vector(int given_layer);
		
		void prepare_sample(bool binary);
		void prepare_sample(bool binary, int given_layer);

		void committ_sample();
	};
	// Comparator
	bool operator <(const Unit& a, const Unit& b);


































	/****************************************
	Visible
	***************************************/

	class UnitVisible : public Unit {

	private:

		// Connections
		std::vector<std::pair<ConnVV*,int>> _conns_vv;
		std::vector<std::pair<ConnVVV*,int>> _conns_vvv;
		std::vector<ConnVH*> _conns_vh;

		// Constructor helpers
		void _clean_up();
		void _move(UnitVisible& other);
		void _copy(const UnitVisible& other);

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
		Sample
		********************/

		void form_propensity_vector();
		void form_propensity_vector(int given_layer);
	};

































	/****************************************
	Hidden
	***************************************/

	class UnitHidden : public Unit {

	private:

		// Layer
		int _layer;

		// Connections
		std::vector<ConnVH*> _conns_vh;
		// Layer -> conns -> (conn,idx of me in conn)
		std::map<int,std::vector<std::pair<ConnHH*,int>>> _conns_hh;

		// Constructor helpers
		void _clean_up();
		void _move(UnitHidden& other);
		void _copy(const UnitHidden& other);

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
		Sample
		********************/

		void form_propensity_vector();
		void form_propensity_vector(int given_layer);
	};
};
