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
     Sampling Options
     ****************************************/
    
    enum class AwakePhaseMode : unsigned int { MEAN_FIELD, GIBBS_SAMPLING };
    
    struct OptionsWakeSleep_BM {
        
        // PCD vs CD
        bool persistent_chains = false;
        
        // Mode to estimate awake phase moments
        AwakePhaseMode awake_phase_mode = AwakePhaseMode::GIBBS_SAMPLING;
        
        // Verbosity
        bool verbose_timing = true;
        
        // Write after awake/asleep
        bool write_after_awake = false;
        bool write_after_asleep = false;
        std::string write_after_awake_dir = "";
        std::string write_after_asleep_dir = "";
    };
    
    struct OptionsWakeSleep_RBM {
      
        // PCD vs CD
        bool persistent_chains = false;
        
        // Verbosity
        bool verbose_timing = true;

    };

	/****************************************
	Lattice
	****************************************/

    typedef std::map<Sptr,arma::vec> layer_occ;
    typedef std::map<int, layer_occ> layers_map;

    class Lattice
	{
    protected:
        
        // ***************
        // MARK: - Misc
        // ***************

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
        // MARK: - Persistent data structures
        // ***************
        
        std::map<int,arma::vec> _pst_prop;
        std::map<int,arma::vec> _pst_r;
        std::map<int,arma::vec> _pst_sign_of_r;
        std::map<int,arma::vec> _pst_sign_of_r_new;
        
        // ***************
        // MARK: Ixn params
        // ***************
        
        // Bias/ixn dicts
        std::vector<Iptr> _all_ixns;
        std::map<int, std::map<Sptr,Iptr>> _bias_dict;
        std::map<int, std::map<Sptr, std::map<int, std::map<Sptr,Iptr>>>> _o2_ixn_dict;

        // ***************
        // MARK: Look up sites
        // ***************
        
        // Lookup a site iterator from x,y,z
        int _look_up_unit(int layer, int x) const;
        int _look_up_unit(int layer, int x, int y) const;
        int _look_up_unit(int layer, int x, int y, int z) const;
        std::vector<int> _look_up_pos(int layer, int idx) const;

	private:

        // ***************
        // MARK: - Misc
        // ***************
        
		// Dimensionality
		int _no_dims;

		// Size
		int _box_length;

        // No markov chains for both awake and asleep
        std::map<MCType,int> _no_markov_chains;
        
        // ***************
        // MARK: - Private methods
        // ***************

        // Binarize
        void _binarize_all_units_in_layer(MCType chain, int i_chain, int layer, bool act);

        // Add ixn to all ixns vec
        void _add_to_all_ixns_vec(Iptr ixn);
        
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
        
        Lattice(int no_dims, int box_length, std::vector<Sptr> species_visible, bool add_visible_layer=true);
        Lattice(const Lattice& other);
		Lattice(Lattice&& other);
		Lattice& operator=(const Lattice& other);
		Lattice& operator=(Lattice&& other);
		virtual ~Lattice();

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
        
        virtual void add_layer(int layer, int box_length, std::vector<Sptr> species);

        // ***************
        // MARK: Biases/ixn params
        // ***************
        
        // Biases
        void set_bias_of_layer(int layer, Sptr sp, Iptr bias);

		// Ixns between layers
        // Always bidirectional
        void set_ixn_between_layers(int layer1, Sptr sp1, int layer2, Sptr sp2, Iptr ixn);

        // Clear biases and ixns
        void clear_all_biases_and_ixns();
        
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
        virtual void set_empty_all_units_in_layer(MCType chain, int i_chain, int layer);
        void set_empty_all_hidden_units(MCType chain, int i_chain);

		// Random
        void set_random_all_units(MCType chain, int i_chain, bool binary);
        void set_random_all_units_in_layer(MCType chain, int i_chain, int layer, bool binary);
        void set_random_all_hidden_units(MCType chain, int i_chain, bool binary);

		// Binarize
        void binarize_all_units(MCType chain, int i_chain);
        virtual void binarize_all_units_in_layer(MCType chain, int i_chain, int layer);
        void binarize_all_hidden_units(MCType chain, int i_chain);

        // ***************
        // MARK: Write/read
        // ***************
        
		virtual void write_layer_to_file(MCType chain, int i_chain, int layer, std::string fname, bool binary) const;
		virtual void read_layer_from_file(MCType chain, int i_chain, int layer, std::string fname, bool binary);
        
        // ***************
        // MARK: Activate layer steps
        // ***************
        
        // For all chains:
        
        // 1. Calculate activation given layer above or below or both
        void activate_layer_calculate_from_below(MCType chain, int layer, double weight_mult=1.0);
        void activate_layer_calculate_from_above(MCType chain, int layer, double weight_mult=1.0);
        void activate_layer_calculate_from_both(MCType chain, int layer, double weight_mult_from_below=1.0, double weight_mult_from_above=1.0);

        // (2) Convert activations to probs
        virtual void activate_layer_convert_to_probs(MCType chain, int layer, bool binary);

        // (3) Commit the new probabilities
        void activate_layer_committ(MCType chain, int layer);

        // ***************
        // MARK: - Reap moments, both awake and asleep
        // ***************
        
        void reap_ixn_moment_diffs() const;
        
        // Manually query moments
        double reap_moment_sample(MCType type, int i_chain, int layer, Sptr species) const;
        double reap_moment_sample(MCType type, int i_chain, int layer_lower, Sptr species_lower, int layer_upper, Sptr species_upper) const;
        double reap_moment(MCType type, int layer, Sptr species) const;
        double reap_moment(MCType type, int layer_lower, Sptr species_lower, int layer_upper, Sptr species_upper) const;
        
        // Query moments for particular ixns
        double reap_moment_sample(MCType type, int i_chain, Iptr ixn) const;
        double reap_moment(MCType type, Iptr ixn) const;

        // Moments for the adjoint terms
        double reap_moment_adjoint_obs_cov_cross_term_sample(int i_chain, Iptr ixn, int layer_domain, Sptr species_domain) const;
        double reap_moment_adjoint_obs_cov_cross_term(Iptr ixn, int layer_domain, Sptr species_domain) const;
                
        // ***************
        // MARK: - Wake/sleep
        // ***************
        
        void wake_sleep_loop_bm(int i_opt_step, int no_steps_awake, int no_steps_asleep, std::vector<FName> &fnames, OptionsWakeSleep_BM options);
        void wake_sleep_loop_rbm(int i_opt_step, int no_cd_steps, std::vector<FName> &fnames, OptionsWakeSleep_RBM options);
        
        // ***************
        // MARK: - Counts
        // ***************
        
        double get_count_1d(MCType chain, int i_chain, Sptr sp) const;
        double get_count_1d(MCType chain, int i_chain, Sptr sp1, Sptr sp2) const;
        double get_count_1d(MCType chain, int i_chain, Sptr sp1, Sptr sp2, Sptr sp3) const;
        double get_count_1d(MCType chain, int i_chain, Sptr sp1, Sptr sp2, Sptr sp3, Sptr sp4) const;
    };
};
