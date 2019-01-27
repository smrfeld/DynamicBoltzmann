#include <string>
#include <vector>
#include <map>

#ifndef FWDS_SPECIES_H
#define FWDS_SPECIES_H
#include "fwds/fwds_species.hpp"
#endif

#include <armadillo>

/************************************
* Namespace for bmla
************************************/

namespace bmla {

	// Forwards
	class UnitVisible;
	class UnitHidden;
	class BiasDict;
	class ConnVH;
	class ConnVV;
	class ConnVVV;
	class ConnHH;
	class O2IxnDict;
	class O3IxnDict;
	class Moment;

	/****************************************
	Lattice
	****************************************/

    typedef std::map<Sptr,arma::vec> layer_occ;
    typedef std::map<int, layer_occ> layers_map;
    
	class Lattice
	{
	private:

		// Dimensionality
		int _no_dims;

		// Size
		int _box_length;

        // No markov chains
        int _no_markov_chains;
        
        // No layers
        int _no_layers;
        
        // Markov chain idx
        // -> Layer idx
        // -> State vector
        std::map<int,layers_map> _latt;
        
		// Hidden layer lookup
		// Layer -> (x,y,z) -> unit 
		std::map<int, std::map<int, int>> _hlookup_1;
		std::map<int, std::map<int, std::map<int, int>>> _hlookup_2;
		std::map<int, std::map<int, std::map<int, std::map<int, int>>>> _hlookup_3;
        
        // Adjacency matrices
        // Idx 0 connects 0, 1
        // Idx i connects i, i+1
        std::map<int,arma::mat> _adj;
        
        // Bias/ixn dicts
        std::map<int,std::shared_ptr<BiasDict>> _bias_dicts;
        std::map<int,std::shared_ptr<O2IxnDict>> _ixn_dicts;
        
		// Species possible
        // layer->species name->species
		std::map<int,std::map<std::string,Sptr>> _species_possible;

		// Lookup a site iterator from x,y,z
		int _look_up_unit(int layer, int x) const;
		int _look_up_unit(int layer, int x, int y) const;
		int _look_up_unit(int layer, int x, int y, int z) const;

        // Add hidden unit to layer
        int _add_hidden_unit(int layer);
        
		// Count helpers
		void _get_count(double &count, Sptr &sp1, Sptr &sp2, const int &idx1, const int &idx2, const int &i_chain, bool binary, bool reversibly) const;
		void _get_count(double &count, Sptr &sp1, Sptr &sp2, Sptr &sp3, const int &idx1, const int &idx2, const int &idx3, const int &i_chain, bool binary, bool reversibly) const;

		// Contructor helpers
		void _clean_up();
		void _move(Lattice& other);
		void _copy(const Lattice& other);

	public:

		/********************
		Constructor
		********************/

        Lattice(int dim, int box_length, std::vector<Sptr> species_visible);
		Lattice(const Lattice& other);
		Lattice(Lattice&& other);
		Lattice& operator=(const Lattice& other);
		Lattice& operator=(Lattice&& other);
		~Lattice();

        /********************
        Getters
         ********************/

        int get_no_dims() const;
        int get_box_length() const;
        
        /********************
        Markov chains
         ********************/

        int get_no_markov_chains() const;
        void set_no_markov_chains(int no_markov_chains);
        
        /********************
        Add a layer
         ********************/

        void add_layer(int layer, int no_units, std::vector<Sptr> species);
        
        /********************
		Biases/ixn params
		********************/

		// Biases
		void set_bias_dict_all_units(std::shared_ptr<BiasDict> bias_dict);
		void set_bias_dict_all_units_in_layer(int layer, std::shared_ptr<BiasDict> bias_dict);

		// Ixns
        void set_ixn_dict_between_layers(int layer_1, int layer_2, std::shared_ptr<O2IxnDict> ixn_dict);

		/********************
		Add connections
		********************/

        void add_conn(int layer1, int x1, int layer2, int x2);
        void add_conn(int layer1, int x1, int y1, int layer2, int x2, int y2);
        void add_conn(int layer1, int x1, int y1, int z1, int layer2, int x2, int y2, int z2);

		/********************
		Add hidden units
		********************/

		void add_hidden_unit(int layer, int x);
		void add_hidden_unit(int layer, int x, int y);
		void add_hidden_unit(int layer, int x, int y, int z);

		/********************
		Apply funcs to all units
		********************/

		// Clear the lattice
        void set_empty_all_units();
        void set_empty_all_units_in_layer(int layer);

		// Random
        void set_random_all_units(bool binary);
        void set_random_all_units_in_layer(int layer, bool binary);

		// Binarize
        void binarize_all_units();
        void binarize_all_units_in_layer(int layer);

		/********************
		Write/read latt to a file
		********************/

		void write_layer_to_file(int layer, std::string fname, bool binary);

		void init_file_reader(std::map<int,std::vector<Sptr>> layers_species_possible);
		void read_layer_from_file(int layer, std::string fname, bool binary);
        
        /********************
        Activate layer
         ********************/
        
        // Activate a specific layer
		void activate_layer(int layer, bool binary);
		void activate_layer(int layer, int given_layer, bool binary);

		/********************
		Get counts for visibles
		********************/

		double get_count_vis(Sptr &sp) const;
		double get_count_vis(Sptr &sp1, Sptr &sp2, bool reversibly) const;
		double get_count_vis(Sptr &sp1, Sptr &sp2, Sptr &sp3, bool reversibly) const;
		double get_count_vis(Sptr &sp1, Sptr &sp2, Sptr &sp3, Sptr &sp4, bool reversibly) const;
	};

};
