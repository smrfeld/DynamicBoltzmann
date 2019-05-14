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
    class Lattice;
    
    enum class MCType: unsigned int;

	/****************************************
	LatticeTraj
	****************************************/

	class LatticeTraj
	{
    protected:
        
        // Lattices
        std::map<int,std::shared_ptr<Lattice>> _lattices;

        // All ixn param trajs
        std::vector<ITptr> _ixn_param_trajs;
        std::map<int, std::map<Sptr,ITptr>> _bias_dict;
        std::map<int, std::map<Sptr, std::map<int, std::map<Sptr,ITptr>>>> _o2_ixn_dict;

    private:
        
        // Add ixn param traj (Unique!)
        void _add_ixn_param_traj(ITptr ixn_param_traj);
        
        // ***************
        // MARK: - Constructor helpers
        // ***************
        
		void _clean_up();
		void _move(LatticeTraj& other);
		void _copy(const LatticeTraj& other);

	public:
        
        // ***************
        // MARK: - Constructor
        // ***************
        
        LatticeTraj(int no_dims, int box_length, std::vector<Sptr> species_visible, bool create_lattice_at_zero=true);
		LatticeTraj(const LatticeTraj& other);
		LatticeTraj(LatticeTraj&& other);
		LatticeTraj& operator=(const LatticeTraj& other);
		LatticeTraj& operator=(LatticeTraj&& other);
		virtual ~LatticeTraj();
        
        // ***************
        // MARK: - No timesteps
        // ***************
        
        // Set no timesteps
        // Idx 0 is timepoint_start_ixn_params in the ixn_params
        virtual void set_no_timesteps(int timepoint_start, int no_timesteps);
        
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
        
        std::shared_ptr<Lattice> get_lattice_at_timepoint(int timepoint) const;
        
        // ***************
        // MARK: - Setup all timepoints
        // ***************
        
        // Add layer
        virtual void add_layer(int layer, int box_length, std::vector<Sptr> species);

		// Biases
        void set_bias_of_layer(int layer, Sptr sp, ITptr bias);

		// Ixns between layers
        // Always bidirectional
        void set_ixn_between_layers(int layer1, Sptr sp1, int layer2, Sptr sp2, ITptr ixn);
        
        // Add connections
        void add_conn(int layer1, int x1, int layer2, int x2);
        void add_conn(int layer1, int x1, int y1, int layer2, int x2, int y2);
        void add_conn(int layer1, int x1, int y1, int z1, int layer2, int x2, int y2, int z2);
    };
};
