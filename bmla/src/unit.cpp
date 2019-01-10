#include "../include/bmla_bits/unit.hpp"

// Other headers
#include "../include/bmla_bits/general.hpp"
#include "../include/bmla_bits/species.hpp"
#include "../include/bmla_bits/ixn_dicts.hpp"
#include "../include/bmla_bits/connections.hpp"

#include <iostream>
#include <fstream>
#include <numeric>
#include "math.h"
#include <ctime>
#include <sstream>
#include <random>
#include <vector>
#include <algorithm>

/************************************
* Namespace for bmla
************************************/

namespace bmla {

	/****************************************
	Class to hold a lattice Unit
	****************************************/

	/********************
	Constructor
	********************/

	Unit::Unit(int x) : Unit(x,0,0) { 
		_dim = 1;
	};
	Unit::Unit(int x, int y) : Unit(x,y,0) {
		_dim = 2;
	};
	Unit::Unit(int x, int y, int z) {
		_dim=3;
		_x = x;
		_y = y;
		_z = z;

		Act act_empty;
		act_empty.sp = nullptr;
		act_empty.act = 0.0;
		act_empty.prob = 1.0;
		act_empty.prop = 1.0;
		_activations.push_back(act_empty);
	};
	Unit::Unit(int x, std::vector<Sptr> species_possible) : Unit(x,0,0,species_possible) { _dim=1; };
	Unit::Unit(int x, int y, std::vector<Sptr> species_possible) : Unit(x,y,0,species_possible) { _dim=2; };
	Unit::Unit(int x, int y, int z, std::vector<Sptr> species_possible) : Unit(x,y,z) {
		set_possible_species(species_possible);
	};	
	Unit::Unit(const Unit& other) {
		_copy(other);
	};
	Unit::Unit(Unit&& other) {
		_move(other);
	};
	Unit& Unit::operator=(const Unit& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	Unit& Unit::operator=(Unit&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
	};
	Unit::~Unit() {
		_clean_up();
	};
	void Unit::_clean_up() {
		// Nothing....
	};
	void Unit::_move(Unit& other) {
		_dim = other._dim;
		_x = other._x;
		_y = other._y;
		_z = other._z;

		_nonzero_occs = other._nonzero_occs;
		_nonzero_occs_tbc = other._nonzero_occs_tbc;

		_sampling_rand = other._sampling_rand;
		_sampling_i_chosen = other._sampling_i_chosen;

		_activations = other._activations;

		_bias_dict = std::move(other._bias_dict);

		other._dim = 0;
		other._x = 0;
		other._y = 0;
		other._z = 0;

		other._nonzero_occs.clear();
		other._nonzero_occs_tbc.clear();

		other._sampling_rand = 0.0;
		other._sampling_i_chosen = 0;

		other._activations.clear();

		other._bias_dict = nullptr;
	};
	void Unit::_copy(const Unit& other) {
		_dim = other._dim;
		_x = other._x;
		_y = other._y;
		_z = other._z;

		_nonzero_occs = other._nonzero_occs;
		_nonzero_occs_tbc = other._nonzero_occs_tbc;

		_sampling_rand = other._sampling_rand;
		_sampling_i_chosen = other._sampling_i_chosen;

		_activations = other._activations;

		_bias_dict = other._bias_dict;
	};

	/********************
	Sample prop vec
	********************/

	void Unit::_sample_prop_vec() {
		// Sample RV
		_sampling_rand = randD(0.0,_activations.back().prop);

		// Empty?
		if (_sampling_rand <= _activations[0].prop) {
			_sampling_i_chosen = 0;
			return;
		};

		// Find interval
		for (int i=1; i<_activations.size(); i++) {
			if (_activations[i-1].prop <= _sampling_rand && _sampling_rand <= _activations[i].prop) {
				_sampling_i_chosen = i;
				return;
			};
		};

		// Never get here
		_sampling_i_chosen = 0;
	};

	/********************
	Print the location and connectivity
	********************/

	void Unit::print() const {
		std::cout << print_str() << std::endl;
	};
	std::string Unit::print_str() const {
		std::ostringstream s;
		if (_dim == 1) {
			s << _x << " ";
		} else if (_dim == 2) {
			s << _x << " " << _y << " ";
		} else if (_dim == 3) {
			s << _x << " " << _y << " " << _z << " ";
		};

		if (_nonzero_occs.size() == 0) {
			s << "empty";
		} else {
			for (auto pr: _nonzero_occs) {
				s << "(" << pr.first->get_name() << "," << pr.second << ") ";
			};
		};
		return s.str();
	};

	/********************
	Location
	********************/

	int Unit::dim() const {
		return _dim;
	};
	int Unit::x() const {
		return _x;
	};
	int Unit::y() const {
		return _y;
	};
	int Unit::z() const {
		return _z;
	};
	bool Unit::less_than(const Unit &other) const {
		if (_dim == 1) {
			return _x < other.x();
		} else if (_dim == 2) {
			return (_x < other.x()) && (_y < other.y());
 		} else if (_dim == 3) {
			return (_x < other.x()) && (_y < other.y()) && (_z < other.z());
	    } else {
	    	return false;
	    };	
	};

	/********************
	Finish setup in lattice
	********************/

	// Add possible species
	void Unit::add_possible_species(Sptr species) {
		// Check that species is not already added
		for (auto const &act: _activations) {
			if (act.sp == species) {
				std::cerr << ">>> Error: Unit::add_possible_species <<< species: " << species->get_name() << " already added." << std::endl;
				exit(EXIT_FAILURE);
			};
		};

		// Add
		Act act;
		act.sp = species;
		act.act = 0.0;
		act.prob = 0.0;
		act.prop = 0.0;
		_activations.push_back(act);
	};
	void Unit::set_possible_species(std::vector<Sptr> species) {
		_nonzero_occs.clear();
		_activations.clear();

		Act act_empty;
		act_empty.sp = nullptr;
		act_empty.act = 0.0;
		act_empty.prob = 1.0;
		act_empty.prop = 1.0;
		_activations.push_back(act_empty);

		for (auto &sp: species) {
			add_possible_species(sp);
		};
	};

	// Bias dict
	void Unit::set_bias_dict(std::shared_ptr<BiasDict> bias_dict) {
		_bias_dict = bias_dict;
	};
	const std::shared_ptr<BiasDict>& Unit::get_bias_dict() const {
		return _bias_dict;
	};

	/********************
	Get probability
	********************/

	double Unit::get_occ(Sptr sp) const {
		if (!sp) {
			double ret = 1.0;
			for (auto pr: _nonzero_occs) {
				ret -= pr.second;
			};
			return ret;
		} else {
			auto it = _nonzero_occs.find(sp);
			if (it != _nonzero_occs.end()) {
				return it->second;
			} else {
				return 0.0;
			};
		};
	};
	const std::unordered_map<Sptr, double>& Unit::get_nonzero_occs() const {
		return _nonzero_occs;
	};
	void Unit::set_occ(Sptr sp, double prob) {
		if (prob > 0.0) {
			_nonzero_occs[sp] = prob;
		};
	};
	void Unit::set_occ_random(bool binary) {
		_nonzero_occs.clear();
		if (binary) {
			int r = randI(1,_activations.size()-1);
			if (r != _activations.size()) {
				_nonzero_occs[_activations[r].sp] = 1.0;
			};
		} else {
			double prob_tot = 0.0;
			for (auto &act: _activations) {
				act.prob = randD(0.0,1.0);
				prob_tot += act.prob;
			};
			for (auto &act: _activations) {
				if (act.sp != nullptr) {
					_nonzero_occs[act.sp] = act.prob/prob_tot;
					// std::cout << "Set: " << act.sp->get_name() << " to: " << _nonzero_occs[act.sp] << std::endl;
				};
			};
		};
	};

	void Unit::binarize() {
		// Check against empty
		if (_nonzero_occs.size() == 0) {
			return; // already binary!
		};

		// Form the propensity vector
		std::vector<double> props;
		std::vector<Sptr> props_sp;
		props.push_back(0.0);
		props_sp.push_back(nullptr);
		// Empty
		double ptot = 0.0;
		for (auto pr: _nonzero_occs) {
			ptot += pr.second;
		};
		props.push_back(1.0-ptot);
		// Others
		for (auto pr: _nonzero_occs) {
			props.push_back(props.back() + pr.second);
			props_sp.push_back(pr.first);
		};

		// Sample RV
		double r = randD(0.0,props.back());

		// Find interval
		for (int i=0; i<props.size()-1; i++) {
			if (props[i] <= r && r <= props[i+1]) {
				// i is the chosen one
				// i = 0 => empty
				// else i is the species chosen in _nonzero_occs
				if (i == 0) {
					set_empty();
				} else {
					set_empty();
					set_occ(props_sp[i],1.0);
				};

				// Done
				return;
			};
		};

		// Never get here
	};

	bool Unit::get_is_empty() const {
		if (_nonzero_occs.size() == 0) {
			return true;
		} else {
			return false;
		};
	};
	void Unit::set_empty() {
		_nonzero_occs.clear();
	};

	/********************
	Get moment
	********************/

	double Unit::get_moment(std::string ixn_param_name) const {
		if (!_bias_dict) {
			std::cerr << ">>> Error: Unit::get_moment <<< no bias dict exists on this unit" << std::endl;
			exit(EXIT_FAILURE);
		};

		double count = 0.0;

		std::vector<Sptr> species = _bias_dict->get_species_from_ixn(ixn_param_name);

		for (auto &sp: species) {
			count += get_occ(sp);
		};

		return count;
	};

	/********************
	Sample
	********************/

	void Unit::form_propensity_vector() {
		for (auto i=0; i<_activations.size(); i++) {
			_activations[i].act = 0.0;
			_activations[i].prob = exp(_activations[i].act);
			_activations[i].prop = _activations[i-1].prop + _activations[i].prob;
		};
	};
	void Unit::form_propensity_vector(int given_layer) {
		for (auto i=0; i<_activations.size(); i++) {
			_activations[i].act = 0.0;
			_activations[i].prob = exp(_activations[i].act);
			_activations[i].prop = _activations[i-1].prop + _activations[i].prob;
		};
	};


	void Unit::prepare_sample(bool binary) {
		// Form the activations vector
		form_propensity_vector();
		// Sample
		_prepare_sample(binary);
	};

	void Unit::prepare_sample(bool binary, int given_layer) {
		// Form the activations vector
		form_propensity_vector(given_layer);
		// Sample
		_prepare_sample(binary);
	};

	void Unit::committ_sample() {
		_nonzero_occs = _nonzero_occs_tbc;
	};

	void Unit::_prepare_sample(bool binary) {

		// Empty the site
		_nonzero_occs_tbc.clear();

		// Commit probs
		if (binary) {

			// Sample RV
			_sample_prop_vec();

			if (_sampling_i_chosen != 0) { // 0 is empty
				// Make the appropriate species at this site
				_nonzero_occs_tbc[_activations[_sampling_i_chosen].sp] = 1.0;
			};

		} else {

			// Write into species
			for (int i=1; i<_activations.size(); i++) {
				_nonzero_occs_tbc[_activations[i].sp] = _activations[i].prob/(_activations.back().prop);
			};
		};
	};
















































	/****************************************
	Class to hold a lattice UnitVisible
	****************************************/

	/********************
	Constructor
	********************/

	UnitVisible::UnitVisible(int x) : Unit(x) {};
	UnitVisible::UnitVisible(int x, int y) : Unit(x,y) {};
	UnitVisible::UnitVisible(int x, int y, int z) : Unit(x,y,z) {};
	UnitVisible::UnitVisible(int x, std::vector<Sptr> species_possible) : Unit(x,species_possible) {};
	UnitVisible::UnitVisible(int x, int y, std::vector<Sptr> species_possible) : Unit(x,y,species_possible) {};
	UnitVisible::UnitVisible(int x, int y, int z, std::vector<Sptr> species_possible) : Unit(x,y,z,species_possible) {};
	UnitVisible::UnitVisible(const UnitVisible& other) : Unit(other) {
		_copy(other);
	};
	UnitVisible::UnitVisible(UnitVisible&& other) : Unit(other) {
		_move(other);
	};
	UnitVisible& UnitVisible::operator=(const UnitVisible& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	UnitVisible& UnitVisible::operator=(UnitVisible&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
	};
	UnitVisible::~UnitVisible() {
		_clean_up();
	};
	void UnitVisible::_clean_up() {
		// Nothing....
	};
	void UnitVisible::_move(UnitVisible& other) {
		_conns_vv = other._conns_vv;
		_conns_vvv = other._conns_vvv;
		_conns_vh = other._conns_vh;

		other._conns_vv.clear();
		other._conns_vvv.clear();
		other._conns_vh.clear();
	};
	void UnitVisible::_copy(const UnitVisible& other) {
		_conns_vv = other._conns_vv;
		_conns_vvv = other._conns_vvv;
		_conns_vh = other._conns_vh;
	};

	/********************
	Add connection
	********************/

	// Add connection
	void UnitVisible::add_conn(ConnVV *conn, int idx_of_me) {
		// Check it does not exist - conns are unique!
		auto pr = std::make_pair(conn,idx_of_me);
		auto it = std::find(_conns_vv.begin(), _conns_vv.end(), pr);
		if (it != _conns_vv.end()) {
			std::cerr << ">>> Error: UnitVisible::add_conn <<< connection already exists!" << std::endl;
			exit(EXIT_FAILURE);
		};
		_conns_vv.push_back(pr);
	};
	const std::vector<std::pair<ConnVV*,int>>& UnitVisible::get_conns_vv() const {
		return _conns_vv;
	};

	void UnitVisible::add_conn(ConnVVV *conn, int idx_of_me) {
		// Check it does not exist - conns are unique!
		auto pr = std::make_pair(conn,idx_of_me);
		auto it = std::find(_conns_vvv.begin(), _conns_vvv.end(), pr);
		if (it != _conns_vvv.end()) {
			std::cerr << ">>> Error: UnitVisible::add_conn <<< connection already exists!" << std::endl;
			exit(EXIT_FAILURE);
		};
		_conns_vvv.push_back(pr);
	};
	const std::vector<std::pair<ConnVVV*,int>>& UnitVisible::get_conns_vvv() const {
		return _conns_vvv;
	};

	void UnitVisible::add_conn(ConnVH *conn) {
		_conns_vh.push_back(conn);
	};
	const std::vector<ConnVH*>& UnitVisible::get_conns_vh() const {
		return _conns_vh;
	};

	/********************
	Sample
	********************/

	void UnitVisible::form_propensity_vector() {
		for (auto i=1; i<_activations.size(); i++) {

			_activations[i].act = 0.0;

			// Bias
			if (_bias_dict) {
				_activations[i].act += _bias_dict->get_ixn(_activations[i].sp);
			};

			// VV conns
			for (auto const &conn_vv: _conns_vv) {
				_activations[i].act += conn_vv.first->get_act_for_species_at_unit(_activations[i].sp,conn_vv.second);
			};

			// VVV conns
			for (auto const &conn_vvv: _conns_vvv) {
				_activations[i].act += conn_vvv.first->get_act_for_species_at_unit(_activations[i].sp,conn_vvv.second);
			};

			// VH conns
			for (auto const &conn_vh: _conns_vh) {
				_activations[i].act += conn_vh->get_act_for_species_at_unit_v(_activations[i].sp);
			};

			_activations[i].prob = exp(_activations[i].act);
			_activations[i].prop = _activations[i-1].prop + _activations[i].prob;
		};
	};
	void UnitVisible::form_propensity_vector(int given_layer) {
		if (given_layer == 1) {

			form_propensity_vector();

		} else {

			// empty
			for (auto i=1; i<_activations.size(); i++) {
				_activations[i].act = 0.0;
				_activations[i].prob = exp(_activations[i].act);
				_activations[i].prop = _activations[i-1].prop + _activations[i].prob;
			};
		};
	};









































































	/****************************************
	Class to hold a lattice UnitHidden
	****************************************/

	/********************
	Constructor
	********************/

	UnitHidden::UnitHidden(int layer, int x) : Unit(x) {
		_layer = layer;
	};
	UnitHidden::UnitHidden(int layer, int x, int y) : Unit(x,y) {
		_layer = layer;
	};
	UnitHidden::UnitHidden(int layer, int x, int y, int z) : Unit(x,y,z) {
		_layer = layer;
	};
	UnitHidden::UnitHidden(const UnitHidden& other) : Unit(other) {
		_copy(other);
	};
	UnitHidden::UnitHidden(UnitHidden&& other) : Unit(other) {
		_move(other);
	};
	UnitHidden& UnitHidden::operator=(const UnitHidden& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	UnitHidden& UnitHidden::operator=(UnitHidden&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
	};
	UnitHidden::~UnitHidden() {
		_clean_up();
	};
	void UnitHidden::_clean_up() {
		// Nothing....
	};
	void UnitHidden::_move(UnitHidden& other) {
		_layer = other._layer;
		_conns_vh = other._conns_vh;
		_conns_hh = other._conns_hh;

		other._layer = 0;
		other._conns_vh.clear();
		other._conns_hh.clear();
	};
	void UnitHidden::_copy(const UnitHidden& other) {
		_layer = other._layer;
		_conns_vh = other._conns_vh;
		_conns_hh = other._conns_hh;
	};

	/********************
	Location
	********************/

	int UnitHidden::layer() const {
		return _layer;
	};

	/********************
	Add connection
	********************/

	// Visible-hidden conns
	void UnitHidden::add_conn(ConnVH *conn) {
		_conns_vh.push_back(conn);
	};
	const std::vector<ConnVH*>& UnitHidden::get_conns_vh() const {
		return _conns_vh;
	};

	// Hidden-hidden conns
	void UnitHidden::add_conn(ConnHH *conn, int idx_of_me, int from_layer) {

		// Pair
		auto pr = std::make_pair(conn,idx_of_me);

		// Find layer
		auto it = _conns_hh.find(from_layer);
		if (it == _conns_hh.end()) {
			// Add
			_conns_hh[from_layer].push_back(pr);
		} else {
			// Check it does not already exist - conns are unique!
			auto it2 = std::find(_conns_hh[from_layer].begin(), _conns_hh[from_layer].end(), pr);
			if (it2 != _conns_hh[from_layer].end()) {
				std::cerr << ">>> Error: UnitHidden::add_conn <<< connection already exists!" << std::endl;
				exit(EXIT_FAILURE);	
			} else {
				// Add
				_conns_hh[from_layer].push_back(pr);
			};
		};
	};
	const std::map<int,std::vector<std::pair<ConnHH*,int>>>& UnitHidden::get_conns_hh() const {
		return _conns_hh;
	};

	/********************
	Sample
	********************/

	void UnitHidden::form_propensity_vector() {

		for (auto i=1; i<_activations.size(); i++) {

			_activations[i].act = 0.0;

			// Bias
			if (_bias_dict) {
				_activations[i].act += _bias_dict->get_ixn(_activations[i].sp);
			};

			// VH conns
			for (auto const &conn_vh: _conns_vh) {
				_activations[i].act += conn_vh->get_act_for_species_at_unit_h(_activations[i].sp);
			};

			// HH conns
			for (auto layer: _conns_hh) {
				// Go through conns
				for (auto const &conn_hh: layer.second) {
					_activations[i].act += conn_hh.first->get_act_for_species_at_unit(_activations[i].sp,conn_hh.second);
				};
			};

			_activations[i].prob = exp(_activations[i].act);
			_activations[i].prop = _activations[i-1].prop + _activations[i].prob;
		};
	};

	void UnitHidden::form_propensity_vector(int given_layer) {

		for (auto i=1; i<_activations.size(); i++) {

			_activations[i].act = 0.0;

			// Bias
			if (_bias_dict) {
				_activations[i].act += _bias_dict->get_ixn(_activations[i].sp);
			};

			// VH conns
			if (given_layer == 0) {
				for (auto const &conn_vh: _conns_vh) {
					_activations[i].act += conn_vh->get_act_for_species_at_unit_h(_activations[i].sp);
				};
			};

			// HH conns
			auto layer = _conns_hh.find(given_layer);
			if (layer != _conns_hh.end()) {
				// Go through conns
				for (auto const &conn_hh: layer->second) {
					_activations[i].act += conn_hh.first->get_act_for_species_at_unit(_activations[i].sp,conn_hh.second);
				};
			};

			_activations[i].prob = exp(_activations[i].act);
			_activations[i].prop = _activations[i-1].prop + _activations[i].prob;
		};	
	};
};
