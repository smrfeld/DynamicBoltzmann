#ifndef IXN_PARAM_TRAJ_h
#define IXN_PARAM_TRAJ_h
#include "ixn_param_traj.hpp"
#endif

#ifndef LATTICE_H
#define LATTICE_H
#include "lattice.hpp"
#endif

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/************************************
	Hidden unit
	************************************/

	// Forward declare
	class Site;

	class HiddenUnit
	{
	private:

		// Number of sites I am connected to
		int _n_conn;

		// Sites I am connected to
		std::vector<Site*> _conn;

		// Ixn param describing the weight connection
		IxnParamTraj *_weight;

		// Time we are at
		int *_t_opt_ptr;

		// Value
		double _val;

		// Activation function
		double _sigma(double x) const;

		// Constructor helpers
		void _clean_up();
		void _copy(const HiddenUnit& other);
		void _copy(HiddenUnit &&other);

	public:

		/********************
		Constructor
		********************/

		HiddenUnit();
		HiddenUnit(std::vector<Site*> conn);
		HiddenUnit(std::vector<Site*> conn, IxnParamTraj *weight);
		HiddenUnit(const HiddenUnit& other);
		HiddenUnit(HiddenUnit&& other);
		HiddenUnit& operator=(const HiddenUnit& other);
		HiddenUnit& operator=(HiddenUnit&& other);
		~HiddenUnit();	

		/********************
		Add connections
		********************/

		void add_conn(Site* site_ptr);

		/********************
		Set pointer to the opt time variable
		********************/

		void set_opt_time_ptr(int *t_opt_ptr);

		/********************
		Activate
		********************/

		void activate(bool binary); 

	};
};