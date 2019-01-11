#include "connection_vh.hpp"

// Other headers
#include "../include/dynamicboltz_bits/general.hpp"
#include "site.hpp"
#include "hidden_unit.hpp"
#include "species.hpp"
#include "ixn_param_traj.hpp"

#include <iostream>
#include <fstream>
#include <numeric>
#include "math.h"
#include <ctime>
#include <sstream>
#include <random>

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	Class to hold a connection from visible to hidden
	****************************************/

	// Constructor
	ConnectionVH::ConnectionVH(Site *site, HiddenUnit *hidden_unit, std::vector<IxnParamTraj*> ips) {
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
	void ConnectionVH::add_ixn_param(IxnParamTraj* ip) {
		// Get the species associated with this ixn param
		std::vector<SpeciesVH> sp_vec = ip->get_species_conn();

		// Go through the species
		for (auto spvh: sp_vec) {
			// Store both ways
			_ips_visible_hidden[spvh.sp_visible][spvh.sp_hidden].push_back(ip);
			_ips_hidden_visible[spvh.sp_hidden][spvh.sp_visible].push_back(ip);
		};
	};

	// Get for a species on the visible unit
	double ConnectionVH::get_act_visible(Species* sp_visible) {
		double act=0.0;
		auto it = _ips_visible_hidden.find(sp_visible);
		if (it != _ips_visible_hidden.end()) {
			// Go through all hidden species
			for (auto iph: it->second) {
				// Go through all ixn params
				for (auto ip: iph.second) {
					// Weight (from ixn param) * hidden units value for this hidden species
					act += ip->get() * _hidden_unit->get_prob(iph.first);
				};
			};
		};
		return act;
	};

	// Get activation for a hidden
	double ConnectionVH::get_act_hidden(HiddenSpecies* sp_hidden) {
		double act=0.0;
		auto it = _ips_hidden_visible.find(sp_hidden);
		if (it != _ips_hidden_visible.end()) {
			// Go through all visible species
			for (auto ipv: it->second) {
				// Go through all ixn params
				for (auto ip: ipv.second) {
					// Weight (from ixn param) * sites value for this visible species
					act += ip->get() * _site->get_prob(ipv.first);
				};
			};
		};
		return act;
	};
};