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

    class LatticeTrajAlternatingBinary : public LatticeTraj
	{
	private:
        
        // Lattices
        std::map<int,std::shared_ptr<LatticeAlternatingBinary>> _lattices_alternating_binary;
        
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
        // MARK: - Get lattice
        // ***************
        
        std::shared_ptr<LatticeAlternatingBinary> get_lattice_alternating_binary_at_timepoint(int timepoint) const;
    };
};
