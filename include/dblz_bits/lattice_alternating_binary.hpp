#include <string>
#include <vector>
#include <map>

#include "fwds/fwds_species.hpp"
#include "fwds/fwds_ixn_param.hpp"

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
	LatticeAlternatingBinary
	****************************************/

    class LatticeAlternatingBinary : public Lattice
	{
	private:

        // Binarize
        void _binarize_all_units_in_layer(MCType chain, int i_chain, int layer, bool act);

        // ***************
        // MARK: - Constructor helpers
        // ***************
        
		void _clean_up();
		void _move(LatticeAlternatingBinary& other);
		void _copy(const LatticeAlternatingBinary& other);

	public:

        // ***************
        // MARK: - Public methods
        // ***************
        
        // ***************
        // MARK: Constructor
        // ***************
        
        LatticeAlternatingBinary(int no_dims, int box_length, std::vector<Sptr> species_visible);
        LatticeAlternatingBinary(const LatticeAlternatingBinary& other);
		LatticeAlternatingBinary(LatticeAlternatingBinary&& other);
		LatticeAlternatingBinary& operator=(const LatticeAlternatingBinary& other);
		LatticeAlternatingBinary& operator=(LatticeAlternatingBinary&& other);
		~LatticeAlternatingBinary();

        // ***************
        // MARK: Clear/Random/Binarize to all units
        // ***************
        
        void set_random_all_units_in_layer(MCType chain, int i_chain, int layer, bool binary);

		// Binarize
        void binarize_all_units_in_layer(MCType chain, int i_chain, int layer);

        // ***************
        // MARK: Write/read
        // ***************
        
		void write_layer_to_file(MCType chain, int i_chain, int layer, std::string fname, bool binary) const;
		void read_layer_from_file(MCType chain, int i_chain, int layer, std::string fname, bool binary);
        
        // ***************
        // MARK: Activate layer steps
        // ***************
        
        // (2) Convert activations to probs
        void activate_layer_convert_to_probs(MCType chain, int layer, bool binary);
    };
};
