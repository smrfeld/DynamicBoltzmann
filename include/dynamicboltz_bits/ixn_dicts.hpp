#ifndef STRING_H
#define STRING_H
#include <string>
#endif

#ifndef MAP_H
#define MAP_H
#include <map>
#endif

#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

#ifndef FWDS_SPECIES_H
#define FWDS_SPECIES_H
#include "fwds/fwds_species.hpp"
#endif

#ifndef FWDS_IXN_PARAM_H
#define FWDS_IXN_PARAM_H
#include "fwds/fwds_ixn_param.hpp"
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	BiasDict
	****************************************/

	class BiasDict {

	private:

		// Site to site to ixns
		std::map<Sptr,std::vector<Iptr>> _dict;

		// Reverse map
		std::map<Iptr,std::vector<Sptr>> _i_dict;
		std::map<std::string,std::vector<Sptr>> _i_str_dict;

		// Constructor helpers
		void _clean_up();
		void _copy(const BiasDict& other);
		void _move(BiasDict& other);

	public:

		// Constructor
		BiasDict();
		BiasDict(const BiasDict& other);
		BiasDict(BiasDict&& other);
		BiasDict& operator=(const BiasDict& other);
		BiasDict& operator=(BiasDict&& other);
		~BiasDict();

		// Add to the dict
		void add_ixn(Sptr sp, Iptr ixn);

		// Get from the dict
		double get_ixn_at_timepoint(Sptr sp, int timepoint) const;

		// Get species
		const std::vector<Sptr>& get_species_from_ixn(Iptr ixn) const;
		const std::vector<Sptr>& get_species_from_ixn(std::string ixn_param_name) const;
	};

	/****************************************
	O2IxnDict
	****************************************/
	
	class O2IxnDict {

	private:

		class Impl;
		std::unique_ptr<Impl> _impl;

	public:

		// Constructor
		O2IxnDict();
		O2IxnDict(const O2IxnDict& other);
		O2IxnDict(O2IxnDict&& other);
		O2IxnDict& operator=(const O2IxnDict& other);
		O2IxnDict& operator=(O2IxnDict&& other);
		~O2IxnDict();

		// Add to the dict
		void add_ixn(Sptr sp_of_site_1, Sptr sp_of_site_2, Iptr ixn, bool reversibly=true);

		// Get from the dict
		double get_ixn_at_timepoint(Sptr sp_of_site_1, Sptr sp_of_site_2, int timepoint) const;

		// Get species
		const std::vector<Sptr2>& get_species_from_ixn(Iptr ixn) const;
		const std::vector<Sptr2>& get_species_from_ixn(std::string ixn_param_name) const;
	};

	/****************************************
	O3IxnDict
	****************************************/
	
	class O3IxnDict {

	private:

		class Impl;
		std::unique_ptr<Impl> _impl;

	public:

		// Constructor
		O3IxnDict();
		O3IxnDict(const O3IxnDict& other);
		O3IxnDict(O3IxnDict&& other);
		O3IxnDict& operator=(const O3IxnDict& other);
		O3IxnDict& operator=(O3IxnDict&& other);
		~O3IxnDict();

		// Add to the dict
		void add_ixn(Sptr sp_of_site_1, Sptr sp_of_site_2, Sptr sp_of_site_3, Iptr ixn, bool reversibly=true);

		// Get from the dict
		double get_ixn_at_timepoint(Sptr sp_of_site_1, Sptr sp_of_site_2, Sptr sp_of_site_3, int timepoint) const;

		// Get species
		const std::vector<Sptr3>& get_species_from_ixn(Iptr ixn) const;
		const std::vector<Sptr3>& get_species_from_ixn(std::string ixn_param_name) const;
	};
};