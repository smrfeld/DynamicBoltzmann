#include <vector>
#include <string>
#include <map>

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
    enum class MCType: unsigned int { AWAKE, ASLEEP };

	/****************************************
	Moment
	****************************************/
	
	class Moment {

	private:

		// Name
		std::string _name;

		// Type = H, J, K, B, W
		IxnParamType _type;

		// No chains
        std::map<MCType, int> _no_markov_chains;

        // Reaped values from the sampler
		// Size = _batch_size
        std::map<MCType, std::vector<double>> _vals_reaped;

		// Averaged values
        std::map<MCType, double> _val_averaged;
	
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

        int get_no_markov_chains(MCType type) const;
        void set_no_markov_chains(MCType type, int no_markov_chains);
        
		/********************
		Reset
		********************/

		void reset_to_zero(MCType type);

		/********************
		Fixed awake
		********************/

		void set_is_awake_moment_fixed(bool flag);
		bool get_is_awake_moment_fixed() const;

		/********************
		Get/set moment
		********************/

		double get_moment(MCType type) const;
		void set_moment(MCType type, double val);

		// Sample
		double get_moment_sample(MCType type, int i_sample) const;
		void set_moment_sample(MCType type, int i_sample, double val);
        void increment_moment_sample(MCType type, int i_sample, double val);

		// Average reaps
		void average_moment_samples(MCType type);

		/********************
		Write
		********************/

		void write_to_file(std::string fname, bool append=false) const;

	};

};
