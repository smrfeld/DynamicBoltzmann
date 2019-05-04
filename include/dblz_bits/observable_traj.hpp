#include <vector>
#include <string>
#include <map>
#include <armadillo>

#include "fwds/fwds_ixn_param_traj.hpp"
#include "fwds/fwds_species.hpp"

/************************************
* Namespace for bmla
************************************/

namespace dblz {

	/****************************************
	ObservableTraj
	****************************************/
	
	class ObservableTraj {

	private:
        
        // Ixn param traj - this is part of the sum
        ITptr _ixn_param;
        
        // The other part of the sum from the domain on the RHS of the diff eq
        int _layer;
        Sptr _species;
        
        // Vals
        int _no_timesteps;
        int _no_timepoints;
        std::vector<double> _vals;

		// Constructor helpers
		void _clean_up();
		void _move(ObservableTraj &other);
		void _copy(const ObservableTraj& other);

	public:

		/********************
		Constructor
		********************/

		ObservableTraj(ITptr ixn_param, int layer, Sptr species);
		ObservableTraj(const ObservableTraj& other);
		ObservableTraj(ObservableTraj&& other);
		ObservableTraj& operator=(const ObservableTraj& other);
		ObservableTraj& operator=(ObservableTraj&& other);
		~ObservableTraj();

        /********************
		Name, type
		********************/

		ITptr get_ixn_param_traj() const;
        int get_layer() const;
        Sptr get_species() const;

        /********************
         Timepoints
         ********************/

        void set_no_timesteps(int no_timesteps);
        int get_no_timesteps() const;
        
		/********************
		Get/set moment
		********************/
        
        // Get moment
        double get_val_at_timepoint(int timepoint) const;
        void increment_val_at_timepoint(int timepoint, double val);
        void set_val_at_timepoint(int timepoint, double val);
        void reset_val_at_timepoint(int timepoint);
	};
};
