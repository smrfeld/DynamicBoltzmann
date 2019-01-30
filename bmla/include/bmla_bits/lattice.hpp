#include <string>
#include <vector>
#include <map>

#ifndef FWDS_SPECIES_H
#define FWDS_SPECIES_H
#include "fwds/fwds_species.hpp"
#endif

#ifndef FWDS_IXN_PARAM_H
#define FWDS_IXN_PARAM_H
#include "fwds/fwds_ixn_param.hpp"
#endif

#include <armadillo>

/************************************
* Namespace for bmla
************************************/

namespace bmla {

	// Forwards
	class Moment;

    enum class MomentType: unsigned int;

	/****************************************
	Lattice
	****************************************/

    typedef std::map<Sptr,arma::vec> layer_occ;
    typedef std::map<int, layer_occ> layers_map;
    typedef std::map<int, std::map<Sptr,std::vector<Iptr>>> bias_dict;
    typedef std::map<int, std::map<Sptr, std::map<int, std::map<Sptr,std::vector<Iptr>>>>> o2_ixn_dict;

	class Lattice
	{
	private:

		// Dimensionality
		int _no_dims;

		// Size
		int _box_length;

        // No markov chains
        // idx 0 = awake stats
        // Other idxs = asleep stats
        int _no_markov_chains;
        int _no_markov_chains_asleep; // = _no_markov_chains - 1

        // No layers
        int _no_layers;
        
        // Current markov chain
        int _i_markov_chain;
        
        // Markov chain idx
        // -> Layer idx
        // -> State vector
        std::map<int,layers_map> _latt;
        std::map<int,layers_map> _latt_act;

		// Layer lookup
		// Layer -> (x,y,z) -> idx
		std::map<int, std::map<int, int>> _lookup_1;
		std::map<int, std::map<int, std::map<int, int>>> _lookup_2;
		std::map<int, std::map<int, std::map<int, std::map<int, int>>>> _lookup_3;
        // Reverse
        // Layer -> idx -> (x,y,z)
        std::map<int, std::map<int, std::vector<int>>> _rlookup;
        
        // Adjacency matrices
        // Idx 0 connects 0, 1
        // Idx i connects i, i+1
        std::map<int,arma::mat> _adj;
        
        // Bias/ixn dicts
        std::vector<Iptr> _all_ixns;
        bias_dict _bias_dict;
        o2_ixn_dict _o2_ixn_dict;
        
        // Multipliers between layers for ixns - NOT bidirectional
        std::map<int, std::map<int, double>> _o2_mults;
        
		// Species possible
        // layer->species name->species
		std::map<int,std::map<std::string,Sptr>> _species_possible;

		// Lookup a site iterator from x,y,z
		int _look_up_unit(int layer, int x) const;
		int _look_up_unit(int layer, int x, int y) const;
		int _look_up_unit(int layer, int x, int y, int z) const;
        std::vector<int> _look_up_pos(int layer, int idx) const;
        
        // Binarize
        void _binarize_all_units_in_layer(int layer, bool act);
        
        // Activations
        void _reset_activations(int layer);
        void _calculate_activations(int layer, int given_layer);
        void _convert_activations(int layer, bool binary);
        
		// Contructor helpers
		void _clean_up();
		void _move(Lattice& other);
		void _copy(const Lattice& other);

	public:

		/********************
		Constructor
		********************/

        Lattice(int no_dims, int box_length, std::vector<Sptr> species_visible);
		Lattice(const Lattice& other);
		Lattice(Lattice&& other);
		Lattice& operator=(const Lattice& other);
		Lattice& operator=(Lattice&& other);
		~Lattice();

        /********************
        Getters
         ********************/

        int get_no_dims() const;
        int get_box_length() const;
        int get_no_units_in_layer(int layer) const;
        
        /********************
        Markov chains
         ********************/

        int get_no_markov_chains_asleep() const;
        void set_no_markov_chains_asleep(int no_markov_chains_asleep);
        void switch_to_markov_chain_asleep(int i_markov_chain_asleep);
        void switch_to_markov_chain_awake();
        
        /********************
        Add a layer
         ********************/

        void add_layer(int layer, int box_length, std::vector<Sptr> species);
        
        /********************
		Biases/ixn params
		********************/

		// Biases
		void add_bias_all_layers(Sptr sp, Iptr bias);
        void add_bias_to_layer(int layer, Sptr sp, Iptr bias);

		// Ixns between layers
        // Always bidirectional
        void add_ixn_between_layers(int layer1, Sptr sp1, int layer2, Sptr sp2, Iptr ixn);

        // Set multiplier
        void set_multiplier_between_layers(int from_layer, int to_layer, double multiplier);
        
        // Get ixns
        double get_bias_in_layer(int layer, Sptr sp) const;
        double get_ixn_between_layers(int from_layer, Sptr from_sp, int to_layer, Sptr to_sp) const;

		/********************
		Add connections
		********************/

        void add_conn(int layer1, int x1, int layer2, int x2);
        void add_conn(int layer1, int x1, int y1, int layer2, int x2, int y2);
        void add_conn(int layer1, int x1, int y1, int z1, int layer2, int x2, int y2, int z2);

		/********************
		Apply funcs to all units
		********************/

		// Clear the lattice
        void set_empty_all_units();
        void set_empty_all_units_in_layer(int layer);
        void set_empty_all_hidden_units();

		// Random
        void set_random_all_units(bool binary);
        void set_random_all_units_in_layer(int layer, bool binary);
        void set_random_all_hidden_units(bool binary);

		// Binarize
        void binarize_all_units();
        void binarize_all_units_in_layer(int layer);
        void binarize_all_hidden_units();

		/********************
		Write/read latt to a file
		********************/

		void write_layer_to_file(int layer, std::string fname, bool binary) const;
		void read_layer_from_file(int layer, std::string fname, bool binary);
        
        /********************
        Activate layer
         ********************/
        
        // Prepare to activate a specific layer
        void activate_layer_prepare(int layer, bool binary);
        void activate_layer_prepare(int layer, int given_layer, bool binary);

        // Commit the activations
        void activate_layer_committ(int layer);

        // Activate a specific layer
        // First prepares
        // Then committs
		void activate_layer(int layer, bool binary);
		void activate_layer(int layer, int given_layer, bool binary);

        // Variational inference
        void variational_inference_hiddens();
        
        // Sample
        void sample(bool binary_visible, bool binary_hidden);
    
        // Make a pass activating upwards
        void activate_upward_pass(bool binary_hidden);
        
		/********************
		Get counts for visibles
		********************/

		double get_count_vis(Sptr &sp) const;
        /*
		double get_count_vis(Sptr &sp1, Sptr &sp2, bool reversibly) const;
		double get_count_vis(Sptr &sp1, Sptr &sp2, Sptr &sp3, bool reversibly) const;
		double get_count_vis(Sptr &sp1, Sptr &sp2, Sptr &sp3, Sptr &sp4, bool reversibly) const;
         */
        
        /********************
        Reap moments
         ********************/

        void reap_moments(MomentType type, int i_sample) const;
	};

};
