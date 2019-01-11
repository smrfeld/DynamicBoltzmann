#include "connection.hpp"

/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Class to hold a connection from visible to hidden
	****************************************/

	// Constructor
	ConnectionVH::ConnectionVH(Site *site, HiddenUnit *hidden_unit, std::vector<IxnParam*> ips) {
		_site = site;
		_hidden_unit = hidden_unit;
		for (auto ip: ips) {
			add_ixn_param(ip);
		};
	};
	ConnectionVH::ConnectionVH(const ConnectionVH& other) {
		_copy(other);
	};
	ConnectionVH::ConnectionVH(ConnectionVH&& other) {
		_copy(other);
		other._reset();
	};
	ConnectionVH& ConnectionVH::operator=(const ConnectionVH& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	ConnectionVH& ConnectionVH::operator=(ConnectionVH&& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
	};
	ConnectionVH::~ConnectionVH() {
		_clean_up();
	};
	void ConnectionVH::_clean_up() {
		// Nothing....
	};
	void ConnectionVH::_reset() {
		_site = nullptr;
		_hidden_unit = nullptr;
		_ips_visible_hidden.clear();
		_ips_hidden_visible.clear();
	};
	void ConnectionVH::_copy(const ConnectionVH& other) {
		_site = other._site;
		_hidden_unit = other._hidden_unit;
		_ips_visible_hidden = other._ips_visible_hidden;
		_ips_hidden_visible = other._ips_hidden_visible;
	};

	// Add ixn param
	void ConnectionVH::add_ixn_param(IxnParam* ip) {
		// Get separately the visible and hidden species associated with this ixn param
		std::vector<Species*> sp_vec_vis = ip->get_species_visible();
		std::vector<Species*> sp_vec_hidd = ip->get_species_hidden();
		// Go through
		for (auto sp_vis: sp_vec_vis) {
			for (auto sp_hidd: sp_vec_hidd) {
				// Add to both
				_ips_visible_hidden[sp_vis][sp_hidd].push_back(ip);
				_ips_hidden_visible[sp_hidd][sp_vis].push_back(ip);
			};
		};
	};

	// Get for a species on the visible unit
	double ConnectionVH::get_act_visible(Species* sp_visible) const {
		double act=0.0;
		auto itv = _ips_visible_hidden.find(sp_visible);
		if (itv != _ips_visible_hidden.end()) {
			// Go through hidden species
			for (auto ith=itv->second.begin(); ith!=itv->second.end(); ith++) {
				// Go through ixn params
				for (auto ip: ith->second) {
					// Weight (from ixn param) * hidden units value
					act += ip->get() * _hidden_unit->get(ith->first);
				};
			};
		};
		return act;
	};

	// Get activation for a hidden
	double ConnectionVH::get_act_hidden(Species *sp_hidden) const {
		double act=0.0;
		auto ith = _ips_hidden_visible.find(sp_hidden);
		if (ith != _ips_hidden_visible.end()) {
			// Go through visible
			for (auto itv = ith->second.begin(); itv != ith->second.end(); itv++) {
				// Go through ixn params
				for (auto ip: itv->second) {
					// Weight (from ixn param) * visible units value
					act += ip->get() * _site->get_prob(itv->first);
				};
			};
		};

		return act;
	};

};