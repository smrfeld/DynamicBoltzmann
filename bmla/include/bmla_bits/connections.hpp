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
	class O2IxnDict;
	class O3IxnDict;

	/****************************************
	ConnVV
	****************************************/

	class ConnVV {
	private:

		UnitVisible *_uv1;
		UnitVisible *_uv2;

		std::shared_ptr<O2IxnDict> _ixn_dict;

		// Constructor helpers
		void _clean_up();
		void _copy(const ConnVV& other);
		void _move(ConnVV& other);

	public:

		/********************
		Constructor
		********************/

		ConnVV(UnitVisible *uv1, UnitVisible *uv2);
		ConnVV(const ConnVV& other);
		ConnVV(ConnVV&& other);
		ConnVV& operator=(const ConnVV& other);
		ConnVV& operator=(ConnVV&& other);
		~ConnVV();

		/********************
		Check units involved
		********************/

		bool check_connects_unit(UnitVisible *uv) const;
		bool check_connects_units(UnitVisible *uv1, UnitVisible *uv2, bool this_order=false) const;

		/********************
		Get count
		********************/

		double get_moment(std::string ixn_param_name, bool binary=true) const;

		/********************
		Set ixns
		********************/

		void set_ixn_dict(std::shared_ptr<O2IxnDict> &ixn_dict);

		/********************
		Get activation on site
		********************/

		// Idx = 1,2
		double get_act_for_species_at_unit(const Sptr &sp_to_place, int idx);
	};

































	
	/****************************************
	ConnVH
	****************************************/

	class ConnVH {
	private:

		UnitVisible *_uv;
		UnitHidden *_uh;

		std::shared_ptr<O2IxnDict> _ixn_dict;

		// Constructor helpers
		void _clean_up();
		void _copy(const ConnVH& other);
		void _move(ConnVH& other);
		
	public:

		/********************
		Constructor
		********************/

		ConnVH(UnitVisible *uv, UnitHidden* uh);
		ConnVH(const ConnVH& other);
		ConnVH(ConnVH&& other);
		ConnVH& operator=(const ConnVH& other);
		ConnVH& operator=(ConnVH&& other);
		~ConnVH();

		/********************
		Check units involved
		********************/

		bool check_connects_unit(UnitVisible *uv) const;
		bool check_connects_unit(UnitHidden *uh) const;
		bool check_connects_units(UnitVisible *uv, UnitHidden *uh) const;

		/********************
		Get count
		********************/

		double get_moment(std::string ixn_param_name, bool binary=true) const;

		/********************
		Set ixns
		********************/

		void set_ixn_dict(std::shared_ptr<O2IxnDict> &ixn_dict);

		/********************
		Get activation on site
		********************/

		double get_act_for_species_at_unit_v(const Sptr &sp_to_place);
		double get_act_for_species_at_unit_h(const Sptr &sp_to_place);
	};
	


































	/****************************************
	ConnVVV
	****************************************/

	class ConnVVV {
	private:

		UnitVisible *_uv1;
		UnitVisible *_uv2;
		UnitVisible *_uv3;

		std::shared_ptr<O3IxnDict> _ixn_dict;

		// Constructor helpers
		void _clean_up();
		void _copy(const ConnVVV& other);
		void _move(ConnVVV& other);

	public:

		/********************
		Constructor
		********************/

		ConnVVV(UnitVisible *uv1, UnitVisible *uv2, UnitVisible *uv3);
		ConnVVV(const ConnVVV& other);
		ConnVVV(ConnVVV&& other);
		ConnVVV& operator=(const ConnVVV& other);
		ConnVVV& operator=(ConnVVV&& other);
		~ConnVVV();

		/********************
		Check units involved
		********************/

		bool check_connects_unit(UnitVisible *uv);
		bool check_connects_units(UnitVisible *uv1, UnitVisible *uv2, bool this_order=false);
		bool check_connects_units(UnitVisible *uv1, UnitVisible *uv2, UnitVisible *uv3, bool this_order=false);

		/********************
		Get moment
		********************/

		double get_moment(std::string ixn_param_name, bool binary=true) const;

		/********************
		Set ixns
		********************/

		void set_ixn_dict(std::shared_ptr<O3IxnDict> &ixn_dict);

		/********************
		Get activation on site
		********************/

		// Idx = 1,2,3
		double get_act_for_species_at_unit(const Sptr &sp_to_place, int idx);
	};
};