#include <string>
#include <vector>
#include <map>

#include "fwds/fwds_species.hpp"
#include "fwds/fwds_ixn_param_traj.hpp"

#include <armadillo>

/************************************
* Namespace for bmla
************************************/

namespace dblz {

	// Forwards
    class LatticeAlternatingBinary;
    
    enum class MCType: unsigned int;

	/****************************************
	LatticeTrajAlternatingBinary
	****************************************/

	class LatticeTrajAlternatingBinary
	{
	private:
        
        // Lattices
        std::map<int,std::shared_ptr<LatticeAlternatingBinary>> _lattices;
        
        // All ixn param trajs
        std::vector<ITptr> _ixn_param_trajs;
        std::map<int, std::map<Sptr,ITptr>> _bias_dict;
        std::map<int, std::map<Sptr, std::map<int, std::map<Sptr,ITptr>>>> _o2_ixn_dict;

        // Add ixn param traj (Unique!)
        void _add_ixn_param_traj(ITptr ixn_param_traj);
        
        // ***************
        // MARK: - Constructor helpers
        // ***************
        
		void _clean_up();
		void _move(LatticeTrajAlternatingBinary& other);
		void _copy(const LatticeTrajAlternatingBinary& other);

	public:
        
        // ***************
        // MARK: - Constructor
        // ***************
        
        LatticeTrajAlternatingBinary(int no_dims, int box_length, std::vector<Sptr> species_visible);
		LatticeTrajAlternatingBinary(const LatticeTrajAlternatingBinary& other);
		LatticeTrajAlternatingBinary(LatticeTrajAlternatingBinary&& other);
		LatticeTrajAlternatingBinary& operator=(const LatticeTrajAlternatingBinary& other);
		LatticeTrajAlternatingBinary& operator=(LatticeTrajAlternatingBinary&& other);
		~LatticeTrajAlternatingBinary();
        
        // ***************
        // MARK: - No timesteps
        // ***************
        
        // Set no timesteps
        // Idx 0 is timepoint_start_ixn_params in the ixn_params
        void set_no_timesteps(int timepoint_start, int no_timesteps);
        
        // ***************
        // MARK: - Getters
        // ***************

        int get_no_dims() const;
        int get_box_length() const;
        int get_no_units_in_layer(int layer) const;
        int get_no_layers() const;
        
        // Get all ixns
        const std::vector<ITptr>& get_all_ixn_param_trajs() const;
        
        // ***************
        // MARK: - Markov chains
        // ***************

        int get_no_markov_chains(MCType type) const;
        void set_no_markov_chains(MCType type, int no_markov_chains);
        
        // ***************
        // MARK: - Get lattice
        // ***************
        
        std::shared_ptr<LatticeAlternatingBinary> get_lattice_at_timepoint(int timepoint) const;
        
        // ***************
        // MARK: - Setup all timepoints
        // ***************
        
        // Add layer
        void add_layer(int layer, int box_length, std::vector<Sptr> species);

		// Biases
        void set_bias_of_layer(int layer, Sptr sp, ITptr bias);

		// Ixns between layers
        // Always bidirectional
        void set_ixn_between_layers(int layer1, Sptr sp1, int layer2, Sptr sp2, ITptr ixn);

        // Set multiplier
        void set_multiplier_between_layers(int from_layer, int to_layer, double multiplier);
        void set_multiplier_for_bias_in_layer(int layer, double multiplier);
        
        // Add connections
        void add_conn(int layer1, int x1, int layer2, int x2);
        void add_conn(int layer1, int x1, int y1, int layer2, int x2, int y2);
        void add_conn(int layer1, int x1, int y1, int z1, int layer2, int x2, int y2, int z2);
    };
};
