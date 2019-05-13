#include <string>
#include <vector>
#include <map>

#include "fwds/fwds_species.hpp"
#include "fwds/fwds_ixn_param_traj.hpp"

#include <armadillo>

#ifndef LATTICE_TRAJ_H
#define LATTICE_TRAJ_H
#include "lattice_traj.hpp"
#endif

#include "fwds/fwds_center_traj.hpp"

/************************************
* Namespace for bmla
************************************/

namespace dblz {

	// Forwards
    class LatticeCenteredHom;
    
    enum class MCType: unsigned int;

	/****************************************
	LatticeTrajCenteredHom
	****************************************/

    class LatticeTrajCenteredHom : public LatticeTraj
	{
	private:
        
        // Lattices
        std::map<int,std::shared_ptr<LatticeCenteredHom>> _lattices_centered_hom;
        
        // Center trajs
        std::map<int,std::vector<CTptr>> _center_trajs;
        
        // ***************
        // MARK: - Constructor helpers
        // ***************
        
		void _clean_up();
		void _move(LatticeTrajCenteredHom& other);
		void _copy(const LatticeTrajCenteredHom& other);

	public:
        
        // ***************
        // MARK: - Constructor
        // ***************
        
        LatticeTrajCenteredHom(int no_dims, int box_length, std::vector<Sptr> species_visible, std::vector<CTptr> center_trajs);
		LatticeTrajCenteredHom(const LatticeTrajCenteredHom& other);
		LatticeTrajCenteredHom(LatticeTrajCenteredHom&& other);
		LatticeTrajCenteredHom& operator=(const LatticeTrajCenteredHom& other);
		LatticeTrajCenteredHom& operator=(LatticeTrajCenteredHom&& other);
		~LatticeTrajCenteredHom();
        
        // ***************
        // MARK: - No timesteps
        // ***************
        
        // Set no timesteps
        // Idx 0 is timepoint_start_ixn_params in the ixn_params
        void set_no_timesteps(int timepoint_start, int no_timesteps);
    
        // ***************
        // MARK: - Setup all timepoints
        // ***************
        
        // Add layer
        void add_layer(int layer, int box_length, std::vector<Sptr> species, std::vector<CTptr> center_trajs);
        // DO NOT USE THE FOLLOWING:
        void add_layer(int layer, int box_length, std::vector<Sptr> species);
        
        // ***************
        // MARK: - Get lattice
        // ***************
        
        std::shared_ptr<LatticeCenteredHom> get_lattice_centered_hom_at_timepoint(int timepoint) const;
    };
};
