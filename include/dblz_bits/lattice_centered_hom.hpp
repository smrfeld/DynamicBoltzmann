#include <string>
#include <vector>
#include <map>

#include "fwds/fwds_species.hpp"
#include "fwds/fwds_ixn_param.hpp"
#include "fwds/fwds_center.hpp"

#include <armadillo>

#ifndef LATTICE_H
#define LATTICE_H
#include "lattice.hpp"
#endif

/************************************
* Namespace for bmla
************************************/

namespace dblz {
    
	/****************************************
	LatticeCenteredHom
	****************************************/

    class LatticeCenteredHom : public Lattice
	{
	private:
        
        // ***************
        // MARK: - Centering
        // ***************
        
        // Connection multiplicities between layers
        // Layer l -> l+1 or l-1 -> conn mult
        std::map<int,std::map<int,int>> _conn_mults;
        
        // Centers
        std::map<int, std::map<Sptr, Cptr>> _centers;
        
        // ***************
        // MARK: - Constructor helpers
        // ***************
        
		void _clean_up();
		void _move(LatticeCenteredHom& other);
		void _copy(const LatticeCenteredHom& other);

	public:

        // ***************
        // MARK: - Public methods
        // ***************
        
        // ***************
        // MARK: Constructor
        // ***************
        
        LatticeCenteredHom(int no_dims, int box_length, std::vector<Sptr> species_visible, std::vector<Cptr> centers);
        LatticeCenteredHom(const LatticeCenteredHom& other);
		LatticeCenteredHom(LatticeCenteredHom&& other);
		LatticeCenteredHom& operator=(const LatticeCenteredHom& other);
		LatticeCenteredHom& operator=(LatticeCenteredHom&& other);
		~LatticeCenteredHom();

        // ***************
        // MARK: Add a layer
        // ***************
        
        // This does nothing
        void add_layer(int layer, int box_length, std::vector<Sptr> species);
        
        // This is the correct function
        void add_layer(int layer, int box_length, std::vector<Sptr> species, std::vector<Cptr> centers, int conn_multiplicity_to_layer_below);

        // ***************
        // MARK: - Centers
        // ***************
        
        void clear_all_centers();
        void set_centers_in_layer(int layer, std::vector<Cptr> centers);
        
        // ***************
        // MARK: - Reap moments, both awake and asleep
        // ***************
        
        void reap_ixn_moment_diffs_and_slide_centers(double sliding_factor, bool calculate_offset=true);

        // ***************
        // MARK: - Write out centers
        // ***************
        
        void read_center_pts_from_file(std::string fname);
        void write_center_pts_to_file(std::string fname) const;
        
        // ***************
        // MARK: - Set centers
        // ***************
        
        Cptr get_center_for_species_in_layer(int layer, Sptr species) const;
    };
};
