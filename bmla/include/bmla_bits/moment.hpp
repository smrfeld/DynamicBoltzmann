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

		// Visible batch size
		int _batch_size;
        
        // Asleep no markov chains
        int _no_markov_chains;

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
		Name, type
		********************/

		std::string get_name() const;
		IxnParamType get_type() const;

		/********************
		Batch size/no markov chains
		********************/

		int get_batch_size() const;
		void set_batch_size(int batch_size);

        int get_no_markov_chains() const;
        void set_no_markov_chains(int no_markov_chains);
        
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

		// Sample
		double get_moment_sample(MomentType type, int i_sample) const;
		void set_moment_sample(MomentType type, int i_sample, double val);
        void increment_moment_sample(MomentType type, int i_sample, double val);

		// Average reaps
		void average_samples(MomentType type);

		/********************
		Write
		********************/

		void write_to_file(std::string fname, bool append=false) const;

	};

};
