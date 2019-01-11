#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

#ifndef MAP_H
#define MAP_H
#include <map>
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	Class to hold a connection from visible to hidden
	****************************************/

	// Forward
	class Site;
	class HiddenUnit;
	class Species;
	class HiddenSpecies;
	class IxnParamTraj;

	class ConnectionVH {
	private:

		// Site
		Site *_site;

		// Hidden unit
		HiddenUnit *_hidden_unit;

		// Ixn param W associated with this connection
		// Stored both directions
		std::map<Species*, std::map<HiddenSpecies*, std::vector<IxnParamTraj*> > > _ips_visible_hidden;
		std::map<HiddenSpecies*, std::map<Species*, std::vector<IxnParamTraj*> > > _ips_hidden_visible;

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const ConnectionVH& other);

	public:

		/********************
		Constructor
		********************/

		ConnectionVH(Site *site, HiddenUnit *hidden_unit, std::vector<IxnParamTraj*> ips);
		ConnectionVH(const ConnectionVH& other);
		ConnectionVH(ConnectionVH&& other);
		ConnectionVH& operator=(const ConnectionVH& other);
		ConnectionVH& operator=(ConnectionVH&& other);
		~ConnectionVH();

		/********************
		Add ixn param
		********************/

		void add_ixn_param(IxnParamTraj* ip);

		/********************
		Get activation for a species on the visible unit
		********************/

		double get_act_visible(Species* sp_visible);

		/********************
		Get activation for a species on the hidden unit
		********************/

		double get_act_hidden(HiddenSpecies* sp_hidden);
	};
};