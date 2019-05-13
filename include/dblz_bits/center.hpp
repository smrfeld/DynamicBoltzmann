#include <string>
#include <vector>
#include "fwds/fwds_species.hpp"

/************************************
 * Namespace for dblz
 ************************************/

namespace dblz {
    
    /****************************************
     Center
     ****************************************/
    
    class Center {
        
    private:
        
        // Value
        double _val;
        double _val_new;
        
        // Layer/species
        int _layer;
        Sptr _species;
        
        // Internal copy func/clean up
        void _clean_up();
        void _copy(const Center& other);
        void _move(Center &other);
        
    public:
        
        // ***************
        // MARK: - Constructor
        // ***************
        
        Center(int layer, Sptr species);
        Center(int layer, Sptr species, double val);
        Center(const Center& other);
        Center& operator=(const Center& other);
        Center(Center&& other);
        Center& operator=(Center&& other);
        ~Center();
        
        // ***************
        // MARK: - Getters/setters
        // ***************
        
        std::string get_name() const;
        
        void set_val(double val);
        double get_val() const;
        
        void reset_val_new();
        void set_val_new(double val_new);
        void increment_val_new(double inc);
        double get_val_new() const;
        
        void slide(double sliding_factor);
        
        // ***************
        // MARK: - Get layer/species
        // ***************
        
        int get_layer() const;
        Sptr get_species() const;
    };
};
