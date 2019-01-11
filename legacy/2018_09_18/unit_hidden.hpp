#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

#ifndef MAP_H
#define MAP_H
#include <map>
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/************************************
	Hidden unit
	************************************/

	// Forward
	class ConnectionVH;
	class HiddenSpecies;
	class IxnParamTraj;

	class UnitHidden
	{
	private:

		// Dimension and location
		int _layer;
		int _dim;
		int _x,_y,_z;

		// Possible species
		std::vector<Sptr> _sp_possible;
		std::unordered_map<std::string, Sptr> _sp_str_map;

		// Connections
		std::vector<ConnectionVH*> _conn;

		// Biases
		std::map<HiddenSpecies*, std::vector<IxnParamTraj*> > _bias;

		// Probs
		std::vector<HiddenSpecies*> _sp_possible;
		std::map<HiddenSpecies*, double> _probs;
		double _prob_empty;

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const UnitHidden& other);

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
		Check setup
		********************/

		void check_setup() const;

		/********************
		Location
		********************/

		int dim() const;
		int layer() const;
		int x() const;
		int y() const;
		int z() const;

		/********************
		Finish setup in lattice
		********************/

		// Add possible species
		void add_possible_species(Sptr species);
		void set_possible_species(std::vector<Sptr> species);
		std::vector<Sptr> get_possible_species() const;

		// Bias dict
		void set_bias_dict(std::shared_ptr<BiasDict> bias_dict);
		std::shared_ptr<BiasDict> get_bias_dict() const;

		// Add a connection
		void add_conn(ConnVH* conn);

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
};