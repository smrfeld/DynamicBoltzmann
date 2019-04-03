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
    class Lattice1DFullyVisible;
    
    enum class MCType: unsigned int;
    enum class LatticeMode: unsigned int;
    
    /****************************************
     LatticeTraj1DFullyVisible
     ****************************************/
    
    class LatticeTraj1DFullyVisible
    {
    private:
        
        // Lattices
        // Time -> Lattice
        std::map<int,std::shared_ptr<Lattice1DFullyVisible>> _lattices;
        
        bool _conn_2, _conn_3;
        
        // All ixn param trajs
        std::vector<ITptr> _ixn_param_trajs;
        std::map<Sptr,ITptr> _bias_dict;
        std::map<Sptr, std::map<Sptr,ITptr>> _o2_ixn_dict;
        std::map<Sptr, std::map<Sptr,std::map<Sptr,ITptr>>> _o3_ixn_dict;
        
        // Add ixn param traj (Unique!)
        void _add_ixn_param_traj(ITptr ixn_param_traj);
        
        // ***************
        // MARK: - Constructor helpers
        // ***************
        
        void _clean_up();
        void _move(LatticeTraj1DFullyVisible& other);
        void _copy(const LatticeTraj1DFullyVisible& other);
        
    public:
        
        // ***************
        // MARK: - Constructor
        // ***************
        
        LatticeTraj1DFullyVisible(int box_length, std::vector<Sptr> species_visible, bool conn_2, bool conn_3);
        LatticeTraj1DFullyVisible(const LatticeTraj1DFullyVisible& other);
        LatticeTraj1DFullyVisible(LatticeTraj1DFullyVisible&& other);
        LatticeTraj1DFullyVisible& operator=(const LatticeTraj1DFullyVisible& other);
        LatticeTraj1DFullyVisible& operator=(LatticeTraj1DFullyVisible&& other);
        ~LatticeTraj1DFullyVisible();
        
        // ***************
        // MARK: - No timesteps
        // ***************
        
        // Set no timesteps
        // Idx 0 is timepoint_start_ixn_params in the ixn_params
        void set_no_timesteps(int timepoint_start, int no_timesteps);
        
        // ***************
        // MARK: - Getters
        // ***************
        
        int get_box_length() const;
        
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
        
        std::shared_ptr<Lattice1DFullyVisible> get_lattice_at_timepoint(int timepoint) const;
        
        // ***************
        // MARK: - Setup all timepoints
        // ***************
        
        // Biases
        void set_bias(Sptr sp, ITptr bias);
        
        // Ixns between layers
        // Always bidirectional
        void set_ixn_2(Sptr sp1, Sptr sp2, ITptr ixn);
        void set_ixn_3(Sptr sp1, Sptr sp2, Sptr sp3, ITptr ixn);
    };
};
