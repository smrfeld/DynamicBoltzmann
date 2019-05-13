#include <string>
#include <vector>
#include "fwds/fwds_species.hpp"
#include "fwds/fwds_center.hpp"

/************************************
 * Namespace for dblz
 ************************************/

namespace dblz {
    
    /****************************************
     CenterTraj
     ****************************************/
    
    class CenterTraj {
        
    private:
        
        // Values
        // Timepoints = timesteps + 1
        int _no_timepoints;
        int _no_timesteps;
        std::vector<Cptr> _centers;
        
        // Layer/species
        int _layer;
        Sptr _species;
        
        // Internal copy func/clean up
        void _clean_up();
        void _copy(const CenterTraj& other);
        void _move(CenterTraj &other);
        
    public:
        
        // ***************
        // MARK: - Constructor
        // ***************
        
        CenterTraj(int layer, Sptr species);
        CenterTraj(const CenterTraj& other);
        CenterTraj& operator=(const CenterTraj& other);
        CenterTraj(CenterTraj&& other);
        CenterTraj& operator=(CenterTraj&& other);
        ~CenterTraj();
        
        // ***************
        // MARK: - Timesteps
        // ***************
        
        void set_no_timesteps(int no_timesteps);
        int get_no_timesteps() const;
        
        // ***************
        // MARK: - Getters/setters
        // ***************
        
        std::string get_name() const;
        
        Cptr get_center_at_timepoint(int timepoint) const;
        
        double get_val_at_timepoint(int timepoint) const;
        
        // Time derivative
        double get_deriv_at_timepoint(int timepoint, double dt, bool bkwd_diff=true) const;
        
        // ***************
        // MARK: - Get layer/species
        // ***************
        
        int get_layer() const;
        Sptr get_species() const;
    };
};
