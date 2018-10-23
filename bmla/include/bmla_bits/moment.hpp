#include <vector>
#include <string>

#ifndef FWDS_SPECIES_H
#define FWDS_SPECIES_H
#include "fwds/fwds_species.hpp"
#endif

/************************************
* Namespace for bmla
************************************/

namespace bmla {

	// Forward
	class UnitVisible;
	class UnitHidden;
	class ConnVV;
	class ConnVVV;
	class ConnVH;
	class ConnHH;
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
		std::vector<ConnHH*> _monitor_x;

		// Batch size
		int _batch_size;

		// Reaped values from the sampler
		// Size = _batch_size
		double *_vals_awake_reaped;
		double *_vals_asleep_reaped;

		// Averaged values
		double _val_awake_averaged;
		double _val_asleep_averaged;
	
		// If the awake moment is fixed
		bool _is_awake_moment_fixed;

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
		void add_conn_to_monitor_x(ConnHH *conn);

		/********************
		Name, type
		********************/

		std::string get_name() const;
		IxnParamType get_type() const;

		/********************
		Batch size
		********************/

		int get_batch_size() const;
		void set_batch_size(int batch_size);

		/********************
		Reset
		********************/

		void reset_to_zero(MomentType type);

		/********************
		Fixed awake
		********************/

		void set_is_awake_moment_fixed(bool flag);
		bool get_is_awake_moment_fixed() const;

		/********************
		Get/set moment
		********************/

		double get_moment(MomentType type) const;
		void set_moment(MomentType type, double val);

		// Batch
		double get_moment_in_batch(MomentType type, int i_batch) const;
		void set_moment_in_batch(MomentType type, int i_batch, double val);

		/********************
		Reap from lattice
		********************/

		void reap_in_batch(MomentType type, int i_batch, bool binary=true);
		
		/********************
		Average reaps
		********************/

		void average_reaps(MomentType type);

		/********************
		Write
		********************/

		void write_to_file(std::string fname, bool append=false) const;

	};

};