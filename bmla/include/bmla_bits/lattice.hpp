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
    typedef std::map<Sptr,std::vector<Iptr>> bias_dict;
    typedef std::map<Sptr,std::map<Sptr,std::vector<Iptr>>> o2_ixn_dict;

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
        std::map<int,bias_dict> _bias_dicts;
        std::map<int,o2_ixn_dict> _o2_ixn_dicts;
        
		// Species possible
        // layer->species name->species
		std::map<int,std::map<std::string,Sptr>> _species_possible;

		// Lookup a site iterator from x,y,z
		int _look_up_unit(int layer, int x) const;
		int _look_up_unit(int layer, int x, int y) const;
		int _look_up_unit(int layer, int x, int y, int z) const;
        std::vector<int> _look_up_pos(int layer, int idx) const;
        
        // Add hidden unit to layer
        int _add_hidden_unit(int layer);
        
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

        Lattice(int dim, int box_length, std::vector<Sptr> species_visible);
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

        // Use idx 0 for awake stats
        // Use other idxs for asleep stats
        int get_no_markov_chains() const;
        void set_no_markov_chains(int no_markov_chains);
        void set_current_markov_chain(int i_markov_chain);
        
        /********************
        Add a layer
         ********************/

        void add_layer(int layer, int no_units, std::vector<Sptr> species);
        
        /********************
		Biases/ixn params
		********************/

		// Biases
		void add_bias_all_layers(Sptr sp, Iptr bias);
        void add_bias_to_layer(int layer, Sptr sp, Iptr bias);

		// Ixns between layer and layer + 1
        void add_ixn_between_layer_and_layer_above(int layer, Sptr sp1, Sptr sp2, Iptr ixn);

        // Get ixns
        double get_bias_in_layer(int layer, Sptr sp) const;
        double get_ixn_between_layer_and_layer_above(int layer, Sptr sp1, Sptr sp2) const;

		/********************
		Add connections
		********************/

        void add_conn(int layer1, int x1, int layer2, int x2);
        void add_conn(int layer1, int x1, int y1, int layer2, int x2, int y2);
        void add_conn(int layer1, int x1, int y1, int z1, int layer2, int x2, int y2, int z2);

		/********************
		Add hidden units
		********************/

		void add_hidden_unit(int layer, int x);
		void add_hidden_unit(int layer, int x, int y);
		void add_hidden_unit(int layer, int x, int y, int z);

		/********************
		Apply funcs to all units
		********************/

		// Clear the lattice
        void set_empty_all_units();
        void set_empty_all_units_in_layer(int layer);

		// Random
        void set_random_all_units(bool binary);
        void set_random_all_units_in_layer(int layer, bool binary);
        void set_random_all_hidden_units(bool binary);

		// Binarize
        void binarize_all_units();
        void binarize_all_units_in_layer(int layer);

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
