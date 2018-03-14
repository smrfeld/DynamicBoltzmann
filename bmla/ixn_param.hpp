#ifndef VECTOR_h
#define VECTOR_h
#include <vector>
#endif

#ifndef STRING_h
#define STRING_h
#include <string>
#endif

#ifndef SPECIES_h
#define SPECIES_h
#include "species.hpp"
#endif

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Interaction parameter
	****************************************/

	// Enumeration of type of dimension
	enum IxnParamType { Hp, Jp };

	class IxnParam {

	private:

		// Name
		std::string _name;

		// Type
		IxnParamType _type;

		// Species
		Species *_sp1;
		Species *_sp2;

		// Number of time points in these trajs
		int _n_t;

		// Value
		double _val;

		// Initial guess
		double _val_guess;

		// Awake and asleep moments
		double _asleep;
		double _awake;

		// Copy, clean up
		void _copy(const IxnParam& other);
		void _copy(IxnParam&& other);
		void _clean_up();

	public:

		// Constructor
		IxnParam(std::string name, IxnParamType type, Species *sp, double val_guess);
		IxnParam(std::string name, IxnParamType type, Species *sp1, Species *sp2, double val_guess);
		IxnParam(const IxnParam& other);
		IxnParam(IxnParam&& other);
		IxnParam & operator=(const IxnParam& other);
		IxnParam & operator=(IxnParam&& other);
		~IxnParam();

		// Update based on diff
		void update(double dopt, bool l2_reg=false, double lambda=0.);

		// Getters/setters
		std::string name() const;
		double get() const;
		void set_guess(double guess);
		void reset();

		// Moments from lattice
		enum MomentType {AWAKE, ASLEEP};
		double get_moment(MomentType moment_type) const;
		void moments_reset(MomentType moment_type);
		void moments_retrieve(MomentType moment_type);
		void moments_retrieve(MomentType moment_type, int batch_size);
		double moments_diff() const;
	};
};

