#include <string>
#include <vector>
#include <map>

#include "fwds/fwds_species.hpp"
#include "fwds/fwds_ixn_param.hpp"

#include <armadillo>

/************************************
* Namespace for bmla
************************************/

namespace dblz {
    
    enum class MCType: unsigned int;
    
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

        // ***************
        // MARK: - Private data
        // ***************
        
		// Dimensionality
		int _no_dims;

		// Size
		int _box_length;

        // No markov chains for both awake and asleep
        std::map<MCType,int> _no_markov_chains;

        // No layers
        int _no_layers;
        
        // ***************
        // MARK: The actual layers
        // ***************
        
        // Markov chain idx
        // -> Layer idx
        // -> State vector
        std::map<MCType,std::map<int,layers_map>> _mc_chains;
        std::map<MCType,std::map<int,layers_map>> _mc_chains_act;
        
        // No units per layer
        std::map<int,int> _no_units_per_layer;
        
		// Layer lookup
		// Layer -> (x,y,z) -> idx
		std::map<int, std::map<int, int>> _lookup_1;
		std::map<int, std::map<int, std::map<int, int>>> _lookup_2;
		std::map<int, std::map<int, std::map<int, std::map<int, int>>>> _lookup_3;
        // Reverse
        // Layer -> idx -> (x,y,z)
        std::map<int, std::map<int, std::vector<int>>> _rlookup;
        
        // Species possible
        // layer->species name->species
        std::map<int,std::map<std::string,Sptr>> _species_possible_map;
        std::map<int,std::vector<Sptr>> _species_possible_vec;

        // ***************
        // MARK: Adjacency matrix
        // ***************
        
        // Adjacency matrices
        // Idx 0 connects 0, 1
        // Idx i connects i, i+1
        std::map<int,arma::mat> _adj;
        
        // ***************
        // MARK: Ixn params
        // ***************
        
        // Bias/ixn dicts
        std::vector<Iptr> _all_ixns;
        bias_dict _bias_dict;
        o2_ixn_dict _o2_ixn_dict;
        
        // Multipliers between layers for ixns - NOT bidirectional
        std::map<int, std::map<int, double>> _o2_mults;
        std::map<int, double> _bias_mults;
        
        // ***************
        // MARK: Batch normalization
        // ***************
        
        // Mode
        bool _bn_mode;
        
        // Batch normalization parameters for the layers
        // (layer) -> (parameter)
        std::map<int,std::map<Sptr,Iptr>> _bn_beta, _bn_gamma;
        
        // Same but for beta-bar and gamma-bar (calculated from beta, gamma, means, vars)
        // (chain type) -> (layer) -> (parameter)
        std::map<MCType,std::map<int,std::map<Sptr,arma::vec>>> _bn_beta_bar, _bn_gamma_bar;
        
        // Means,vars over the current batch
        // CHANGES over the steps of mean field / gibbs sampling
        // (chain type) -> (layer) -> means/vars
        std::map<MCType,std::map<int,std::map<Sptr,arma::vec>>> _bn_means, _bn_vars;
        
        // Small epsilon to prevent div by zero
        double _bn_eps;
        
        // ***************
        // MARK: - Private methods
        // ***************

        // Binarize
        void _binarize_all_units_in_layer(MCType chain, int i_chain, int layer, bool act);

        // ***************
        // MARK: Look up sites
        // ***************
        
		// Lookup a site iterator from x,y,z
		int _look_up_unit(int layer, int x) const;
		int _look_up_unit(int layer, int x, int y) const;
		int _look_up_unit(int layer, int x, int y, int z) const;
        std::vector<int> _look_up_pos(int layer, int idx) const;
        
        // ***************
        // MARK: Batch normalization
        // ***************
        
        // Calculate means from activations
        void _bn_calculate_means_vars_from_activations(MCType chain, int layer);
        
        // Calculate bar parameters from means, vars
        void _bn_calculate_bar_params_from_means_vars(MCType chain, int layer);
        
        // Apply BN transform
        void _bn_apply_affine_transform_to_all_chains(MCType chain, int layer);
        
        // ***************
        // MARK: Activations
        // ***************
        
        // Reset activations to 0
        void _reset_activations(MCType chain, int i_chain, int layer);
        
        // Calculate bias
        void _calculate_bias(MCType chain, int i_chain, int layer);
        
        // Calculate activation given layer above or below
        void _calculate_activations_from_below(MCType chain, int i_chain, int layer);
        void _calculate_activations_from_above(MCType chain, int i_chain, int layer);

        // Calculate activations when the scale paramaters gamma in eqn (7) are not incorporated into th weights
        // See bullet pt 2 on p. 366
        void _calculate_activations_from_above_bn(MCType chain, int i_chain, int layer);
        void _calculate_activations_from_below_bn(MCType chain, int i_chain, int layer);

        // ***************
        // MARK: - Constructor helpers
        // ***************
        
		void _clean_up();
		void _move(Lattice& other);
		void _copy(const Lattice& other);

	public:

        // ***************
        // MARK: - Public methods
        // ***************
        
        // ***************
        // MARK: Constructor
        // ***************
        
        Lattice(int no_dims, int box_length, std::vector<Sptr> species_visible, bool batch_norm_mode);
		Lattice(const Lattice& other);
		Lattice(Lattice&& other);
		Lattice& operator=(const Lattice& other);
		Lattice& operator=(Lattice&& other);
		~Lattice();

        // ***************
        // MARK: Getters
        // ***************

        int get_no_dims() const;
        int get_box_length() const;
        int get_no_units_in_layer(int layer) const;
        int get_no_layers() const;
        
        // ***************
        // MARK: Markov chains
        // ***************

        int get_no_markov_chains(MCType type) const;

        void set_no_markov_chains(MCType type, int no_markov_chains);
        
        // ***************
        // MARK: Add a layer
        // ***************
        
        void add_layer(int layer, int box_length, std::vector<Sptr> species);
        void add_layer(int layer, int box_length, std::vector<Sptr> species, Iptr beta, Iptr gamma);

        // ***************
        // MARK: Biases/ixn params
        // ***************
        
		// Biases
		void add_bias_all_layers(Sptr sp, Iptr bias);
        void add_bias_to_layer(int layer, Sptr sp, Iptr bias);

		// Ixns between layers
        // Always bidirectional
        void add_ixn_between_layers(int layer1, Sptr sp1, int layer2, Sptr sp2, Iptr ixn);

        // Set multiplier
        void set_multiplier_between_layers(int from_layer, int to_layer, double multiplier);
        void set_multiplier_for_bias_in_layer(int layer, double multiplier);

        // Get ixns
        double get_bias_in_layer(int layer, Sptr sp) const;
        double get_ixn_between_layers(int from_layer, Sptr from_sp, int to_layer, Sptr to_sp) const;

        // Get all ixns
        const std::vector<Iptr>& get_all_ixn_params() const;
        
        // ***************
        // MARK: Add connections
        // ***************

        void add_conn(int layer1, int x1, int layer2, int x2);
        void add_conn(int layer1, int x1, int y1, int layer2, int x2, int y2);
        void add_conn(int layer1, int x1, int y1, int z1, int layer2, int x2, int y2, int z2);

        // ***************
        // MARK: Clear/Random/Binarize to all units
        // ***************

		// Clear the lattice
        void set_empty_all_units(MCType chain, int i_chain);
        void set_empty_all_units_in_layer(MCType chain, int i_chain, int layer);
        void set_empty_all_hidden_units(MCType chain, int i_chain);

		// Random
        void set_random_all_units(MCType chain, int i_chain, bool binary);
        void set_random_all_units_in_layer(MCType chain, int i_chain, int layer, bool binary);
        void set_random_all_hidden_units(MCType chain, int i_chain, bool binary);

		// Binarize
        void binarize_all_units(MCType chain, int i_chain);
        void binarize_all_units_in_layer(MCType chain, int i_chain, int layer);
        void binarize_all_hidden_units(MCType chain, int i_chain);

        // ***************
        // MARK: Write/read
        // ***************
        
		void write_layer_to_file(MCType chain, int i_chain, int layer, std::string fname, bool binary) const;
		void read_layer_from_file(MCType chain, int i_chain, int layer, std::string fname, bool binary);
        
        // ***************
        // MARK: Activate layer steps
        // ***************
        
        // For all chains:
        
        // (1.a) Calculate activations for a specific layer
        // Both directions
        void activate_layer_calculate(MCType chain, int layer);
        // Only one direction
        void activate_layer_calculate(MCType chain, int layer, int given_layer);
        // (1.b) Alternatively, include batch norm params
        // NOTE: For these two, BN params must already exist! See below to calculate!
        void activate_layer_calculate_bn(MCType chain, int layer);
        void activate_layer_calculate_bn(MCType chain, int layer, int given_layer);
        
        // (2) Convert activations to probs
        void activate_layer_convert_to_probs(MCType chain, int layer, bool binary);

        // (3) Commit the new probabilities
        void activate_layer_committ(MCType chain, int layer);

        // ***************
        // MARK: - Calculate params in BN mode
        // ***************
        
        // Calculate BN params
        void calculate_bn_params(MCType chain);
        
        // ***************
        // MARK: - Mean field / gibbs sampling
        // ***************
        
        // Variational inference ie mean field
        void mean_field_hiddens_step();

        // Gibbs sampling
        void gibbs_sampling_step(bool binary_visible, bool binary_hidden);
        void gibbs_sampling_step_parallel(bool binary_visible, bool binary_hidden);

        // Make a pass activating upwards
        void activate_upward_pass(MCType chain, bool binary_hidden);
        void activate_upward_pass_with_2x_weights_1x_bias(MCType chain, bool binary_hidden);

        // ***************
        // MARK: Get counts for visible layer
        // ***************

        /*
        double get_count_vis(Sptr &sp) const;
		double get_count_vis(Sptr &sp1, Sptr &sp2, bool reversibly) const;
		double get_count_vis(Sptr &sp1, Sptr &sp2, Sptr &sp3, bool reversibly) const;
		double get_count_vis(Sptr &sp1, Sptr &sp2, Sptr &sp3, Sptr &sp4, bool reversibly) const;
         */
        
        // ***************
        // MARK: Reap moments
        // ***************

        void reap_moments(MCType chain);
	};
};
