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
    struct FName;
    
    /****************************************
     OptProblemStatic Options
     ****************************************/
    
    struct OptionsWakeSleep {
        
        // Verbosity
        bool verbose = false;
        bool verbose_timing = true;
        
        // Sampling options
        // Is the visible reconstruction binary, EXCEPT in the last phase?
        bool is_asleep_visible_binary = true;
        // Is the hidden reconstruction binary, EXCEPT in the last phase?
        bool is_asleep_hidden_binary = true;
        // Is the visible reconstruction binary in the last phase?
        bool is_asleep_visible_binary_final = true;
        // Is the hidden reconstruction binary in the last phase?
        bool is_asleep_hidden_binary_final = false;
        
        // Write after awake/asleep
        bool write_after_awake = false;
        bool write_after_asleep = false;
        std::string write_after_awake_dir = "";
        std::string write_after_asleep_dir = "";
        
        // Gibbs sampling awake phase
        bool gibbs_sample_awake_phase = false;
        bool gibbs_sample_awake_phase_hidden_binary = true;
        
        // Sliding factor
        double sliding_factor = 0.01;
    };
    
	/****************************************
	Lattice
	****************************************/

    typedef std::map<Sptr,arma::vec> layer_occ;
    typedef std::map<int, layer_occ> layers_map;

    // First option
    enum class LatticeMode: unsigned int { NORMAL, NORMAL_W_CENTERED_GRADIENT_VEC, NORMAL_W_CENTERED_GRADIENT_PT, CENTERED_PT };

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
        // Layer 1 -> layer 2 -> matrix
        // corresponds to:
        // layer 2 . matrix . layer 1
        // i.e. "layer 1 TO layer 2" 
        std::map<int,std::map<int,arma::sp_mat>> _adj;
        
        // ***************
        // MARK: Ixn params
        // ***************
        
        // Bias/ixn dicts
        std::vector<Iptr> _all_ixns;
        std::map<int, std::map<Sptr,Iptr>> _bias_dict;
        std::map<int, std::map<Sptr, std::map<int, std::map<Sptr,Iptr>>>> _o2_ixn_dict;
        
        // Multipliers between layers for ixns - NOT bidirectional
        std::map<int, std::map<int, double>> _o2_mults;
        std::map<int, double> _bias_mults;
        
        // ***************
        // MARK: Mode: normal vs centered vs batch normalized
        // ***************
        
        // Mode
        LatticeMode _mode;
        
        // ***************
        // MARK: - Centering
        // ***************
        
        // Moving average means for each layer
        std::map<int, std::map<Sptr, arma::vec>> _c_means;
        std::map<int, std::map<Sptr, arma::vec>> _c_batch_means;
        
        // As above, but pointwise
        std::map<int, std::map<Sptr, double>> _cpt_means;
        std::map<int, std::map<Sptr, double>> _cpt_batch_means;
        
        // ***************
        // MARK: - Private methods
        // ***************

        // Binarize
        void _binarize_all_units_in_layer(MCType chain, int i_chain, int layer, bool act);

        // Add ixn to all ixns vec
        void _add_to_all_ixns_vec(Iptr ixn);
        
        // ***************
        // MARK: Look up sites
        // ***************
        
		// Lookup a site iterator from x,y,z
		int _look_up_unit(int layer, int x) const;
		int _look_up_unit(int layer, int x, int y) const;
		int _look_up_unit(int layer, int x, int y, int z) const;
        std::vector<int> _look_up_pos(int layer, int idx) const;
        
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
        
        Lattice(int no_dims, int box_length, std::vector<Sptr> species_visible, LatticeMode mode);
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
        LatticeMode get_lattice_mode() const;
        
        // ***************
        // MARK: Markov chains
        // ***************

        int get_no_markov_chains(MCType type) const;
        void set_no_markov_chains(MCType type, int no_markov_chains);
        
        // ***************
        // MARK: Add a layer
        // ***************
        
        void add_layer(int layer, int box_length, std::vector<Sptr> species);

        // ***************
        // MARK: Biases/ixn params
        // ***************
        
        // Biases
        void set_bias_of_layer(int layer, Sptr sp, Iptr bias);

		// Ixns between layers
        // Always bidirectional
        void set_ixn_between_layers(int layer1, Sptr sp1, int layer2, Sptr sp2, Iptr ixn);

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
        
        // 1. Calculate activation given layer above or below or both
        void activate_layer_calculate_from_below(MCType chain, int layer); // figures out which normal or centered
        void activate_layer_calculate_from_below_centered(MCType chain, int layer);
        void activate_layer_calculate_from_below_normal(MCType chain, int layer);

        void activate_layer_calculate_from_above(MCType chain, int layer); // figures out which normal or centered
        void activate_layer_calculate_from_above_centered(MCType chain, int layer);
        void activate_layer_calculate_from_above_normal(MCType chain, int layer);

        void activate_layer_calculate_from_both(MCType chain, int layer); // figures out which normal or centered
        void activate_layer_calculate_from_both_centered(MCType chain, int layer);
        void activate_layer_calculate_from_both_normal(MCType chain, int layer);

        // (2) Convert activations to probs
        void activate_layer_convert_to_probs(MCType chain, int layer, bool binary);

        // (3) Commit the new probabilities
        void activate_layer_committ(MCType chain, int layer);

        // ***************
        // MARK: - Mean field / gibbs sampling
        // ***************
        
        // Variational inference ie mean field
        void mean_field_hiddens_step(); // figures out which normal or centered
        void mean_field_hiddens_step_normal();
        void mean_field_hiddens_step_centered();

        // Gibbs sampling for awake phase
        void gibbs_sampling_step_awake(bool binary_hidden); // figures out which normal or centered
        void gibbs_sampling_step_awake_normal(bool binary_hidden);
        void gibbs_sampling_step_awake_centered(bool binary_hidden);

        void gibbs_sampling_step_parallel_awake(bool binary_hidden); // figures out which normal or centered
        void gibbs_sampling_step_parallel_awake_normal(bool binary_hidden);
        void gibbs_sampling_step_parallel_awake_centered(bool binary_hidden);

        // Gibbs sampling
        void gibbs_sampling_step(bool binary_visible, bool binary_hidden); // figures out which normal or centered
        void gibbs_sampling_step_normal(bool binary_visible, bool binary_hidden);
        void gibbs_sampling_step_centered(bool binary_visible, bool binary_hidden);

        void gibbs_sampling_step_parallel(bool binary_visible, bool binary_hidden); // figures out which normal or centered
        void gibbs_sampling_step_parallel_normal(bool binary_visible, bool binary_hidden);
        void gibbs_sampling_step_parallel_centered(bool binary_visible, bool binary_hidden);

        // Make a pass activating upwards
        void activate_upward_pass(MCType chain, bool binary_hidden); // figures out which normal or centered
        void activate_upward_pass_normal(MCType chain, bool binary_hidden);
        void activate_upward_pass_centered(MCType chain, bool binary_hidden);

        void activate_upward_pass_with_2x_weights_1x_bias(MCType chain, bool binary_hidden); // figures out which normal or centered

        // ***************
        // MARK: - Reap moments, both awake and asleep
        // ***************
        
        void reap_moments_and_slide_centers_normal();
        void reap_moments_and_slide_centers_normal_w_centered_gradient_vec(bool slide_means, double sliding_factor);
        void reap_moments_and_slide_centers_normal_w_centered_gradient_pt(bool slide_means, double sliding_factor);
        void reap_moments_and_slide_centers_centered_pt(bool slide_means, double sliding_factor);

        // ***************
        // MARK: - Write out centers
        // ***************
        
        void read_centers_from_file(int layer, std::string fname);
        void write_centers_to_file(int layer, std::string fname) const;
        
        void read_center_pts_from_file(std::string fname);
        void write_center_pts_to_file(std::string fname) const;
        
        // ***************
        // MARK: - Set centers
        // ***************
        
        double get_center_for_species_in_layer(int layer, Sptr species) const;
        void set_center_for_species_in_layer(int layer, Sptr species, double center);

        // ***************
        // MARK: - Wake/sleep
        // ***************
        
        void wake_sleep_loop(int i_opt_step, int no_mean_field_updates, int no_gibbs_sampling_steps, std::vector<FName> &fnames, OptionsWakeSleep options);
	};
};
