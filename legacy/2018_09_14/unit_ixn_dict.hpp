#ifndef MAP_H
#define MAP_H
#include <map>
#endif

#ifndef FWDS_SPECIES_H
#define FWDS_SPECIES_H
#include "fwds_species.hpp"
#endif

#ifndef FWDS_IXN_PARAM_H
#define FWDS_IXN_PARAM_H
#include "fwds_ixn_param.hpp"
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	UnitIxnDict
	****************************************/

	class UnitIxnDict {

	private:

		// Site to site to ixns
		std::map<Sptr,std::vector<Iptr>> _dict;

		// Constructor helpers
		void _clean_up();
		void _copy(const UnitIxnDict& other);
		void _move(UnitIxnDict& other);

	public:

		// Constructor
		UnitIxnDict();
		UnitIxnDict(const UnitIxnDict& other);
		UnitIxnDict(UnitIxnDict&& other);
		UnitIxnDict& operator=(const UnitIxnDict& other);
		UnitIxnDict& operator=(UnitIxnDict&& other);
		~UnitIxnDict();

		// Add to the dict
		void add_ixn(Sptr sp, Iptr ixn);

		// Get from the dict
		double get_ixn(Sptr sp, int timepoint) const;
	};

};