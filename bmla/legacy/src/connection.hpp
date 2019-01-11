/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Class to hold a connection from visible to hidden
	****************************************/

	class ConnectionVH {
	private:

		// Site
		Site *_site;

		// Hidden unit
		HiddenUnit *_hidden_unit;

		// Ixn params
		// Store both by: vis species -> hidden species
		// and: hidden species -> vis species
		std::map<Species*, std::map<Species*, std::vector<IxnParam*> > > _ips_visible_hidden;
		std::map<Species*, std::map<Species*, std::vector<IxnParam*> > > _ips_hidden_visible;

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const ConnectionVH& other);

	public:

		// Constructor
		ConnectionVH(Site *site, HiddenUnit *hidden_unit, std::vector<IxnParam*> ips);
		ConnectionVH(const ConnectionVH& other);
		ConnectionVH(ConnectionVH&& other);
		ConnectionVH& operator=(const ConnectionVH& other);
		ConnectionVH& operator=(ConnectionVH&& other);
		~ConnectionVH();

		// Add ixn param
		void add_ixn_param(IxnParam* ip);

		// Get for a species on the visible unit
		double get_act_visible(Species* sp_visible) const;

		// Get activation for a hidden 
		double get_act_hidden(Species *sp_hidden) const;
	};
};