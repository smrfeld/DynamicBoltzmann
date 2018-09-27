#include <vector>
#include <string>

#ifndef FWDS_SPECIES_H
#define FWDS_SPECIES_H
#include "fwds/fwds_species.hpp"
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	// Forward
	class UnitVisible;
	class UnitHidden;
	class ConnVV;
	class ConnVVV;
	class ConnVH;
	enum class IxnParamType: unsigned int;

	/****************************************
	Moment
	****************************************/
	
	enum class MomentType: unsigned int { AWAKE, ASLEEP };

	class Moment {

	private:

		// Name
		std::string _name;

		// Type = H, J, K, B, W
		IxnParamType _type;

		// Units to monitor
		std::vector<UnitVisible*> _monitor_h;
		std::vector<UnitHidden*> _monitor_b;
		std::vector<ConnVV*> _monitor_j;
		std::vector<ConnVVV*> _monitor_k;
		std::vector<ConnVH*> _monitor_w;

		// No time points = time steps + 1
		int _no_timesteps;
		int _no_timepoints;

		// Batch size
		int _batch_size;

		// Reaped values from the sampler
		// Size = _no_timesteps * _batch_size
		double *_vals_awake_reaped;
		double *_vals_asleep_reaped;

		// Averaged values
		// Size = _no_timesteps
		double *_vals_awake_averaged;
		double *_vals_asleep_averaged;
	
		// Constructor helpers
		void _clean_up();
		void _move(Moment &other);
		void _copy(const Moment& other);

	public:

		/********************
		Constructor
		********************/

		Moment(std::string name, IxnParamType type);
		Moment(const Moment& other);
		Moment(Moment&& other);
		Moment& operator=(const Moment& other);
		Moment& operator=(Moment&& other);
		~Moment();

		/********************
		Verbose
		********************/

		void print_moment_comparison() const;

		/********************
		Finish setup
		********************/

		// Add units
		void add_unit_to_monitor_h(UnitVisible *uv);
		void add_unit_to_monitor_b(UnitHidden *uh);
		void add_conn_to_monitor_j(ConnVV *conn);
		void add_conn_to_monitor_k(ConnVVV *conn);
		void add_conn_to_monitor_w(ConnVH *conn);

		/********************
		Name, type
		********************/

		std::string get_name() const;
		IxnParamType get_type() const;

		/********************
		Number timesteps
		********************/

		int get_no_timesteps() const;
		void set_no_timesteps(int no_timesteps);

		/********************
		Batch size
		********************/

		int get_batch_size() const;
		void set_batch_size(int batch_size);

		/********************
		Reset
		********************/

		void reset_to_zero();

		/********************
		Get/set moment
		********************/

		double get_moment_at_timepoint(MomentType type, int timepoint) const;
		void set_moment_at_timepoint(MomentType type, int timepoint, double val);

		// Batch
		double get_moment_at_timepoint_in_batch(MomentType type, int timepoint, int i_batch) const;
		void set_moment_at_timepoint_in_batch(MomentType type, int timepoint, int i_batch, double val);

		/********************
		Reap from lattice
		********************/

		void reap_as_timepoint_in_batch(MomentType type, int timepoint, int i_batch, bool binary=true);
		
		/********************
		Average reaps
		********************/

		void average_reaps_as_timepoint(MomentType type, int timepoint);

		/********************
		Write
		********************/

		void write_to_file(std::string fname) const;

	};

};