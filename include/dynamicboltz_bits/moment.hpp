#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

#ifndef STRING_H
#define STRING_H
#include <string>
#endif

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
	class ConnVV;
	class ConnVVV;
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
		std::vector<ConnVV*> _monitor_j;
		std::vector<ConnVVV*> _monitor_k;

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

		// Species
		/*
		std::vector<Sptr> _sp_h;
		std::vector<Sptr> _sp_b;
		std::vector<Sptr2> _sp_j;
		std::vector<Sptr3> _sp_k;
		std::vector<Sptr2> _sp_w;
		*/
	
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
		Check setup
		********************/

		void check_setup() const;

		/********************
		Verbose
		********************/

		void print_moment_comparison() const;

		/********************
		Finish setup
		********************/

		// Add units
		void add_unit_to_monitor_h(UnitVisible *uv);
		void add_conn_to_monitor_j(ConnVV *conn);
		void add_conn_to_monitor_k(ConnVVV *conn);

		// Add species
		/*
		void add_species_h(Sptr species);
		void add_species_b(Sptr species);
		void add_species_j(Sptr species1, Sptr species2);
		void add_species_k(Sptr species1, Sptr species2, Sptr species3);
		void add_species_w(Sptr speciesV, Sptr speciesH);
		*/

		/********************
		Get species
		********************/

		/*
		const std::vector<Sptr>& get_species_h() const;
		const std::vector<Sptr>& get_species_b() const;
		const std::vector<Sptr2>& get_species_j() const;
		const std::vector<Sptr3>& get_species_k() const;
		const std::vector<Sptr2>& get_species_w() const;
		*/

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
	};

};