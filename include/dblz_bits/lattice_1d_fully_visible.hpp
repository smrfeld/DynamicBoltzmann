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
    
    struct OptionsWakeSleep_1DFV_CD {
      
        // Verbosity
        bool verbose_timing = true;

    };

	/****************************************
	Lattice
	****************************************/

    typedef std::map<Sptr,arma::vec> layer_occ;
    typedef std::map<int, layer_occ> layers_map;

    class Lattice1DFullyVisible
	{
	private:

        // ***************
        // MARK: - Private data
        // ***************
        
		// Size
		int _box_length;

        // No markov chains for both awake and asleep
        std::map<MCType,int> _no_markov_chains;

        // ***************
        // MARK: The actual layers
        // ***************
        
        // Markov chain idx
        // -> State vector
        std::map<MCType,layers_map> _mc_chains;
        std::map<MCType,layers_map> _mc_chains_act;
        
		// Layer lookup
		// x -> idx
		std::map<int, int> _lookup_1;
        // Reverse
        // idx -> x
        std::map<int, int> _rlookup;
        
        // Species possible
        // species name->species
        std::map<std::string,Sptr> _species_possible_map;
        std::vector<Sptr> _species_possible_vec;

        // ***************
        // MARK: Connectivity
        // ***************
        
        bool _conn_2; // nns
        bool _conn_3; // triplets
        
        // ***************
        // MARK: Ixn params
        // ***************
        
        // Bias/ixn dicts
        std::vector<Iptr> _all_ixns;
        std::map<Sptr,Iptr> _bias_dict;
        std::map<Sptr,std::map<Sptr,Iptr>> _o2_ixn_dict;
        std::map<Sptr,std::map<Sptr,std::map<Sptr,Iptr>>> _o3_ixn_dict;
        
        // ***************
        // MARK: - Private methods
        // ***************

        // Add ixn to all ixns vec
        void _add_to_all_ixns_vec(Iptr ixn);
        
        // ***************
        // MARK: - Persistent data structures
        // ***************
        
        arma::vec _pst_prop;
        arma::vec _pst_r;
        arma::vec _pst_sign_of_r;
        arma::vec _pst_sign_of_r_new;

        // ***************
        // MARK: Look up sites
        // ***************
        
		// Lookup a site iterator from x,y,z
		int _look_up_unit(int x) const;
        int _look_up_pos(int idx) const;
        
        // ***************
        // MARK: - Constructor helpers
        // ***************
        
		void _clean_up();
		void _move(Lattice1DFullyVisible& other);
		void _copy(const Lattice1DFullyVisible& other);

	public:

        // ***************
        // MARK: - Public methods
        // ***************
        
        // ***************
        // MARK: Constructor
        // ***************
        
        Lattice1DFullyVisible(int box_length, std::vector<Sptr> species_visible, bool conn_2, bool conn_3);
        Lattice1DFullyVisible(const Lattice1DFullyVisible& other);
		Lattice1DFullyVisible(Lattice1DFullyVisible&& other);
		Lattice1DFullyVisible& operator=(const Lattice1DFullyVisible& other);
		Lattice1DFullyVisible& operator=(Lattice1DFullyVisible&& other);
		~Lattice1DFullyVisible();

        // ***************
        // MARK: Getters
        // ***************

        int get_box_length() const;
        int get_no_units() const;
        
        // ***************
        // MARK: Markov chains
        // ***************

        int get_no_markov_chains(MCType type) const;
        void set_no_markov_chains(MCType type, int no_markov_chains);
        
        // ***************
        // MARK: Biases/ixn params
        // ***************
        
        // Biases
        void set_bias(Sptr sp, Iptr bias);

        // Always bidirectional
        void set_ixn_2(Sptr sp1, Sptr sp2, Iptr ixn);
        void set_ixn_3(Sptr sp1, Sptr sp2, Sptr sp3, Iptr ixn);

        // Get ixns
        double get_bias(Sptr sp) const;
        double get_ixn_2(Sptr sp1, Sptr sp2) const;
        double get_ixn_3(Sptr sp1, Sptr sp2, Sptr sp3) const;

        // Get all ixns
        const std::vector<Iptr>& get_all_ixn_params() const;
        
        // Clear biases and ixns
        void clear_all_biases_and_ixns();

        // ***************
        // MARK: Clear/Random/Binarize to all units
        // ***************

		// Clear the lattice
        void set_empty_all_units(MCType chain, int i_chain);

		// Random
        void set_random_all_units(MCType chain, int i_chain, bool binary);

		// Binarize
        void binarize_unit(MCType chain, int i_chain, int x, bool act);
        void binarize_all_units(MCType chain, int i_chain, bool act);

        // ***************
        // MARK: Write/read
        // ***************
        
		void write_to_file(MCType chain, int i_chain, std::string fname, bool binary) const;
		void read_from_file(MCType chain, int i_chain, std::string fname, bool binary);
        
        // ***************
        // MARK: Activate layer steps
        // ***************
        
        // For all chains:
        
        // 1. Calculate activation given layer above or below or both
        void activate_calculate(MCType chain, int idx_start, int skip);

        // (2) Convert activations to probs
        void activate_convert_to_probs(MCType chain, bool binary, int idx_start, int skip);

        // (3) Commit the new probabilities
        void activate_committ(MCType chain, int idx_start, int skip);

        // ***************
        // MARK: - Mean field / gibbs sampling
        // ***************
        
        // Gibbs sampling
        void gibbs_sampling_step();

        // ***************
        // MARK: - Reap moments, both awake and asleep
        // ***************
        
        void reap_ixn_moment_diffs();

        // ***************
        // MARK: - Wake/sleep
        // ***************
        
        void wake_sleep_loop_cd(int i_opt_step, int no_cd_steps, std::vector<FName> &fnames, OptionsWakeSleep_1DFV_CD options);
        
        // ***************
        // MARK: - Counts
        // ***************
        
        double get_count(MCType chain, int i_chain, Sptr sp) const;
        double get_count(MCType chain, int i_chain, Sptr sp1, Sptr sp2) const;
        double get_count(MCType chain, int i_chain, Sptr sp1, Sptr sp2, Sptr sp3) const;
        double get_count(MCType chain, int i_chain, Sptr sp1, Sptr sp2, Sptr sp3, Sptr sp4) const;
    };
};
