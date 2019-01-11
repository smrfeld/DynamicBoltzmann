#include "../../include/dynamicboltz_bits/unit_visible.hpp"

// Other headers
#include "../../include/dynamicboltz_bits/general.hpp"
#include "../../include/dynamicboltz_bits/species.hpp"
#include "../../include/dynamicboltz_bits/ixn_dicts.hpp"
#include "../../include/dynamicboltz_bits/connections.hpp"

#include <iostream>
#include <fstream>
#include <numeric>
#include "math.h"
#include <ctime>
#include <sstream>
#include <random>
#include <vector>

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	Class to hold a lattice site
	***************************************/

	class UnitVisible::Impl {

	private:

		// Dimensionality and location
		int _dim;
		int _x,_y,_z;

		// Possible species
		std::vector<Sptr> _sp_possible;
		std::unordered_map<std::string, Sptr> _sp_str_map;

		// Binary mode
		bool _b_mode_flag;
		Sptr _b_mode_sp; // occupying species

		// Probabilistic mode
		std::unordered_map<std::string, double*> _p_mode_probs_str;
		std::unordered_map<Sptr, double> _p_mode_probs;
		double _p_mode_prob_empty;

		// Bias dict
		std::shared_ptr<BiasDict> _bias_dict;

		// Connections
		std::vector<std::pair<ConnVV*,int>> _conns_vv;
		std::vector<std::pair<ConnVVV*,int>> _conns_vvv;
		std::vector<ConnVH*> _conns_vh;

		// Data structures for sampling
		double _sampling_prob, _sampling_rand, _sampling_tot;
		int _sampling_i_chosen;
		std::vector<double> _sampling_probs, _sampling_props;

		// Sample a vector of propensities (cumulative probabilities)
		void _sample_prop_vec();

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const Impl& other);

	public:

		/********************
		Constructor
		********************/

		Impl(int x);
		Impl(int x, std::vector<Sptr> species_possible);
		Impl(int x, int y);
		Impl(int x, int y, std::vector<Sptr> species_possible);
		Impl(int x, int y, int z);
		Impl(int x, int y, int z, std::vector<Sptr> species_possible);
		Impl(const Impl& other);
		Impl(Impl&& other);
		Impl& operator=(const Impl& other);
		Impl& operator=(Impl&& other);
		~Impl();

		/********************
		Check setup
		********************/

		void check_setup() const;

		/********************
		Location
		********************/

		int dim() const;
		int x() const;
		int y() const;
		int z() const;
		bool less_than(const UnitVisible &other) const;

		/********************
		Finish setup in lattice
		********************/

		// Add possible species
		void add_possible_species(Sptr species);
		void set_possible_species(std::vector<Sptr> species);
		const std::vector<Sptr>& get_possible_species() const;

		// Bias dict
		void set_bias_dict(std::shared_ptr<BiasDict> bias_dict);
		const std::shared_ptr<BiasDict>& get_bias_dict() const;

		// Add connection
		void add_conn(ConnVV *conn, int idx_of_me);
		std::vector<std::pair<ConnVV*,int>> get_conns_vv() const;

		void add_conn(ConnVVV *conn, int idx_of_me);
		std::vector<std::pair<ConnVVV*,int>> get_conns_vvv() const;
		
		void add_conn(ConnVH *conn);

		/********************
		Get probability
		********************/

		// Check mode
		bool check_is_b_mode() const;

		// Flip between the two modes
		void set_b_mode(bool flag);

		// Binary
		Sptr get_b_mode_species() const; // nullptr for empty
		void set_b_mode_species(Sptr sp);
		void set_b_mode_species(std::string sp);
		void set_b_mode_empty();
		bool check_b_mode_is_empty() const;

		// Probabilistic
		double get_p_mode_prob(Sptr sp) const; // nullptr for empty
		const std::unordered_map<Sptr, double>& get_p_mode_probs() const;
		void set_p_mode_prob(Sptr sp, double prob);
		void set_p_mode_prob(std::string sp, double prob);
		void set_p_mode_empty();
		bool check_p_mode_is_empty() const;

		// Convert between the two
		void convert_b_to_p_mode();
		void convert_p_to_b_mode();

		/********************
		Get moment
		********************/

		double get_moment(std::string ixn_param_name, bool binary=true) const;

		/********************
		Get activation
		********************/

		double get_act_for_species_at_timepoint(Sptr &species, int timepoint) const;

		/********************
		Sample
		********************/

		void sample_at_timepoint(int timepoint, bool binary);
	};









































	/****************************************
	Class to hold a lattice UnitVisible
	****************************************/

	/********************
	Constructor
	********************/

	UnitVisible::Impl::Impl(int x) : Impl(x,0,0) { 
		_dim = 1;
	};
	UnitVisible::Impl::Impl(int x, int y) : Impl(x,y,0) {
		_dim = 2;
	};
	UnitVisible::Impl::Impl(int x, int y, int z) {
		_dim=3;
		_x = x;
		_y = y;
		_z = z;
		_p_mode_prob_empty = 1.0; // default = empty

		// Default = binary
		_b_mode_flag = true;
		_b_mode_sp = nullptr;

		// Setup for sampling
		_sampling_props.push_back(0.0); // 1st element always 0 for probs
		_sampling_props.push_back(1.0); // empty species
		_sampling_probs.push_back(1.0); // empty species
	};
	UnitVisible::Impl::Impl(int x, std::vector<Sptr> species_possible) : Impl(x,0,0,species_possible) { _dim=1; };
	UnitVisible::Impl::Impl(int x, int y, std::vector<Sptr> species_possible) : Impl(x,y,0,species_possible) { _dim=2; };
	UnitVisible::Impl::Impl(int x, int y, int z, std::vector<Sptr> species_possible) : Impl(x,y,z) {
		set_possible_species(species_possible);
	};	
	UnitVisible::Impl::Impl(const Impl& other) {
		_copy(other);
	};
	UnitVisible::Impl::Impl(Impl&& other) {
		_copy(other);
		other._reset();
	};
	UnitVisible::Impl& UnitVisible::Impl::operator=(const Impl& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	UnitVisible::Impl& UnitVisible::Impl::operator=(Impl&& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
	};
	UnitVisible::Impl::~Impl() {
		_clean_up();
	};
	void UnitVisible::Impl::_clean_up() {
		// Nothing....
	};
	void UnitVisible::Impl::_reset() {
		_dim = 0;
		_x = 0;
		_y = 0;
		_z = 0;

		_sp_possible.clear();
		_sp_str_map.clear();

		_b_mode_flag = true;
		_b_mode_sp = nullptr;

		_p_mode_probs_str.clear();
		_p_mode_probs.clear();
		_p_mode_prob_empty = 0.0;

		_bias_dict = nullptr;

		_conns_vv.clear();
		_conns_vvv.clear();
		_conns_vh.clear();

		_sampling_prob = 0.0;
		_sampling_rand = 0.0;
		_sampling_tot = 0.0;
		_sampling_i_chosen = 0;
		_sampling_props.clear();
		_sampling_probs.clear();
	};
	void UnitVisible::Impl::_copy(const Impl& other) {
		_dim = other._dim;
		_x = other._x;
		_y = other._y;
		_z = other._z;

		_sp_possible = other._sp_possible;
		_sp_str_map = other._sp_str_map;

		_b_mode_flag = other._b_mode_flag;
		_b_mode_sp = other._b_mode_sp;

		_p_mode_probs_str = other._p_mode_probs_str;
		_p_mode_probs = other._p_mode_probs;
		_p_mode_prob_empty = other._p_mode_prob_empty;

		_bias_dict = other._bias_dict;

		_conns_vv = other._conns_vv;
		_conns_vvv = other._conns_vvv;
		_conns_vh = other._conns_vh;

		_sampling_prob = other._sampling_prob;
		_sampling_rand = other._sampling_rand;
		_sampling_tot = other._sampling_tot;
		_sampling_i_chosen = other._sampling_i_chosen;
		_sampling_props = other._sampling_props;
		_sampling_probs = other._sampling_probs;
	};

	/********************
	Sample prop vec
	********************/

	void UnitVisible::Impl::_sample_prop_vec() {
		// Sample RV
		_sampling_rand = randD(0.0,_sampling_props.back());

		// Find interval
		for (int i=0; i<_sampling_props.size()-1; i++) {
			if (_sampling_props[i] <= _sampling_rand && _sampling_rand <= _sampling_props[i+1]) {
				_sampling_i_chosen = i;
				return;
			};
		};
		_sampling_i_chosen = 0; // never get here
	};

	/********************
	Print the location and connectivity
	********************/

	void UnitVisible::Impl::check_setup() const {
		// std::cout << _x << " " << _y << " " << _z << " hidden nbrs: (" << _hidden_conns.size() << ") visible nbrs: (" << _nbrs.size() << ")" << std::endl;
	};

	/********************
	Location
	********************/

	int UnitVisible::Impl::dim() const {
		return _dim;
	};
	int UnitVisible::Impl::x() const {
		return _x;
	};
	int UnitVisible::Impl::y() const {
		return _y;
	};
	int UnitVisible::Impl::z() const {
		return _z;
	};
	bool UnitVisible::Impl::less_than(const UnitVisible &other) const {
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
	void UnitVisible::Impl::add_possible_species(Sptr species) {
		// Check that species is not already added
		auto it = std::find(_sp_possible.begin(), _sp_possible.end(), species);
		if (it != _sp_possible.end()) {
			std::cerr << ">>> Error: UnitVisible::Impl::add_possible_species <<< species: " << species->get_name() << " already added." << std::endl;
			exit(EXIT_FAILURE);
		};

		// Add
		_sp_possible.push_back(species);
		_sp_str_map[species->get_name()] = species;

		// Probs
		_p_mode_probs[species] = 0.0;
		_p_mode_probs_str[species->get_name()] = &(_p_mode_probs[species]);

		// Sampling
		_sampling_probs.push_back(0.0);
		_sampling_props.push_back(0.0);
	};
	void UnitVisible::Impl::set_possible_species(std::vector<Sptr> species) {
		_p_mode_probs.clear();
		_sp_possible.clear();
		_sp_str_map.clear();

		_sampling_probs.clear();
		_sampling_props.clear();

		// Empty for sampling
		_sampling_props.push_back(0.0); // 1st element always 0 for probs
		_sampling_props.push_back(1.0); // empty species
		_sampling_probs.push_back(1.0); // empty species

		for (auto &sp: species) {
			add_possible_species(sp);
		};
	};
	const std::vector<Sptr>& UnitVisible::Impl::get_possible_species() const {
		return _sp_possible;
	};

	// Bias dict
	void UnitVisible::Impl::set_bias_dict(std::shared_ptr<BiasDict> bias_dict) {
		_bias_dict = bias_dict;
	};
	const std::shared_ptr<BiasDict>& UnitVisible::Impl::get_bias_dict() const {
		return _bias_dict;
	};

	// Add connection
	void UnitVisible::Impl::add_conn(ConnVV *conn, int idx_of_me) {
		// Check it does not exist - conns are unique!
		auto pr = std::make_pair(conn,idx_of_me);
		auto it = std::find(_conns_vv.begin(), _conns_vv.end(), pr);
		if (it != _conns_vv.end()) {
			std::cerr << ">>> Error: UnitVisible::Impl::add_conn <<< connection already exists!" << std::endl;
			exit(EXIT_FAILURE);
		};
		_conns_vv.push_back(pr);
	};
	std::vector<std::pair<ConnVV*,int>> UnitVisible::Impl::get_conns_vv() const {
		return _conns_vv;
	};

	void UnitVisible::Impl::add_conn(ConnVVV *conn, int idx_of_me) {
		// Check it does not exist - conns are unique!
		auto pr = std::make_pair(conn,idx_of_me);
		auto it = std::find(_conns_vvv.begin(), _conns_vvv.end(), pr);
		if (it != _conns_vvv.end()) {
			std::cerr << ">>> Error: UnitVisible::Impl::add_conn <<< connection already exists!" << std::endl;
			exit(EXIT_FAILURE);
		};
		_conns_vvv.push_back(pr);
	};
	std::vector<std::pair<ConnVVV*,int>> UnitVisible::Impl::get_conns_vvv() const {
		return _conns_vvv;
	};

	void UnitVisible::Impl::add_conn(ConnVH *conn) {
		_conns_vh.push_back(conn);
	};

	/********************
	Get probability
	********************/

	// Check mode
	bool UnitVisible::Impl::check_is_b_mode() const {
		if (_b_mode_flag) {
			return true;
		} else {
			return false;
		};
	};

	// Flip between the two modes
	void UnitVisible::Impl::set_b_mode(bool flag) {
		_b_mode_flag = flag;
	};

	// Binary
	Sptr UnitVisible::Impl::get_b_mode_species() const {
		return _b_mode_sp;
	};
	void UnitVisible::Impl::set_b_mode_species(Sptr sp) {
		_b_mode_sp = sp;
	};
	void UnitVisible::Impl::set_b_mode_species(std::string sp) {
		_b_mode_sp = _sp_str_map[sp];
	};
	void UnitVisible::Impl::set_b_mode_empty() {
		_b_mode_sp = nullptr;
	};
	bool UnitVisible::Impl::check_b_mode_is_empty() const {
		if (_b_mode_sp) {
			return false;
		} else {
			return true;
		};
	};

	// Probabilistic
	double UnitVisible::Impl::get_p_mode_prob(Sptr sp) const {
		// nullptr for prob of empty
		if (sp == nullptr) {
			return _p_mode_prob_empty;
		};

		auto it = _p_mode_probs.find(sp);
		if (it != _p_mode_probs.end()) {
			return it->second;
		} else {
			return 0.0;
		};
	};
	const std::unordered_map<Sptr, double>& UnitVisible::Impl::get_p_mode_probs() const {
		return _p_mode_probs;
	};
	void UnitVisible::Impl::set_p_mode_prob(Sptr sp, double prob) {
		// nullptr for prob of empty
		if (sp == nullptr) {
			_p_mode_prob_empty = prob;
			return;
		};

		_p_mode_probs[sp] = prob;
	};
	void UnitVisible::Impl::set_p_mode_prob(std::string sp, double prob) {
		*(_p_mode_probs_str[sp]) = prob;
	};
	void UnitVisible::Impl::set_p_mode_empty() {
		// Remove counts on existing species
		for (auto pr: _p_mode_probs) {
			if (pr.second > 0.0) {
				// Clear
				_p_mode_probs[pr.first] = 0.0;
			};
		};

		// Empty prob = 1
		_p_mode_prob_empty = 1.0;
	};
	bool UnitVisible::Impl::check_p_mode_is_empty() const {
		if (_p_mode_prob_empty == 1.0) {
			return true;
		} else {
			return false;
		};
	};

	// Convert between the two
	void UnitVisible::Impl::convert_b_to_p_mode() {
		// Clear current
		set_p_mode_empty();

		// Occupancy
		if (_b_mode_sp) {
			set_p_mode_prob(_b_mode_sp,1.0);
		};

		// Binary mode = false
		_b_mode_flag = false;
	};
	void UnitVisible::Impl::convert_p_to_b_mode() {

		// 0th element is always 0
		// _sampling_props.push_back(0.0);

		// 1st element=different; reset at end to 1
		_sampling_props[1] = _p_mode_prob_empty;

		// Others
		for (auto i=0; i<_sp_possible.size(); i++) {
			_sampling_props[i+2] = _sampling_props[i+1] + _p_mode_probs[_sp_possible[i]];
		};

		// Sample
		_sample_prop_vec();

		if (_sampling_i_chosen==0) {
			// Flip down (new spin = 0)
			set_b_mode_empty();
		} else {
			// Make the appropriate species at this site (guaranteed empty)
			set_b_mode_species(_sp_possible[_sampling_i_chosen-1]);
		};

		// Reset propensity of empty to 1
		_sampling_props[1] = 1.0;

		// Binary mode = true
		_b_mode_flag = true;
	};

	/********************
	Get moment
	********************/

	double UnitVisible::Impl::get_moment(std::string ixn_param_name, bool binary) const {
		if (!_bias_dict) {
			std::cerr << ">>> Error: UnitVisible::Impl::get_moment <<< no bias dict exists on this unit" << std::endl;
			exit(EXIT_FAILURE);
		};

		double count = 0.0;

		std::vector<Sptr> species = _bias_dict->get_species_from_ixn(ixn_param_name);
		for (auto &sp: species) {

			if (binary) {

				// Binary
				if (_b_mode_sp == sp) {
					count += 1.0;
				};

			} else {

				// Prob
				count += get_p_mode_prob(sp);
			};
		};

		return count;
	};

	/********************
	Get activation
	********************/

	double UnitVisible::Impl::get_act_for_species_at_timepoint(Sptr &species, int timepoint) const {
		double act = 0.0;

		// Bias
		if (_bias_dict) {
			act += _bias_dict->get_ixn_at_timepoint(species,timepoint);
		};

		// VV conns
		for (auto const &conn_vv: _conns_vv) {
			act += conn_vv.first->get_act_for_species_at_unit_at_timepoint(species,conn_vv.second,timepoint);
		};

		// VVV conns
		for (auto const &conn_vvv: _conns_vvv) {
			act += conn_vvv.first->get_act_for_species_at_unit_at_timepoint(species,conn_vvv.second,timepoint);
		};

		return act;
	};

	/********************
	Sample
	********************/

	void UnitVisible::Impl::sample_at_timepoint(int timepoint, bool binary) {
		// 0 element of props is always 0
		// _sampling_props.push_back(0.0);

		// Empty = 1
		// 1st element of props = 0th element of probs is always 1
		// _sampling_props.push_back(1.0);
		// _sampling_probs.push_back(1.0);

		// Commit probs
		if (binary) {

			// Go through all possible species this could be, calculate propensities
			for (auto i=0; i<_sp_possible.size(); i++) {
				
				_sampling_prob = exp(get_act_for_species_at_timepoint(_sp_possible[i],timepoint));

				// Append prop
				_sampling_props[i+2] = _sampling_props[i+1] + _sampling_prob;
			};

			// Sample RV
			_sample_prop_vec();

			if (_sampling_i_chosen==0) {
				// Flip down (new spin = 0)
				set_b_mode_empty();
			} else {
				// Make the appropriate species at this site (guaranteed empty)
				set_b_mode_species(_sp_possible[_sampling_i_chosen-1]);
			};

		} else {

			// Go through all possible species this could be, calculate propensities
			_sampling_tot = 0.0;
			for (auto i=0; i<_sp_possible.size(); i++) {
				
				_sampling_prob = exp(get_act_for_species_at_timepoint(_sp_possible[i],timepoint));

				// Append prob
				_sampling_probs[i+1] = _sampling_prob;
				_sampling_tot += _sampling_prob;
			};

			// Write into species
			Sptr s_empty = nullptr;
			set_p_mode_prob(s_empty,_sampling_probs[0]/_sampling_tot);
			for (int i=0; i<_sp_possible.size(); i++) {
				set_p_mode_prob(_sp_possible[i],_sampling_probs[i+1]/_sampling_tot);
			};
		};
	};















































	/****************************************
	Forwards
	****************************************/

	/********************
	Constructor
	********************/

	UnitVisible::UnitVisible(int x) : _impl(new Impl(x)) {};
	UnitVisible::UnitVisible(int x, std::vector<Sptr> species_possible) : _impl(new Impl(x,species_possible)) {};
	UnitVisible::UnitVisible(int x, int y) : _impl(new Impl(x,y)) {};
	UnitVisible::UnitVisible(int x, int y, std::vector<Sptr> species_possible) : _impl(new Impl(x,y,species_possible)) {};
	UnitVisible::UnitVisible(int x, int y, int z) : _impl(new Impl(x,y,z)) {};
	UnitVisible::UnitVisible(int x, int y, int z, std::vector<Sptr> species_possible) : _impl(new Impl(x,y,z,species_possible)) {};
	UnitVisible::UnitVisible(const UnitVisible& other) : _impl(new Impl(*other._impl)) {};
	UnitVisible::UnitVisible(UnitVisible&& other) : _impl(std::move(other._impl)) {};
	UnitVisible& UnitVisible::operator=(const UnitVisible& other) {
        _impl.reset( new Impl( *other._impl ) );
        return *this; 
	};
	UnitVisible& UnitVisible::operator=(UnitVisible&& other) {
        _impl = std::move(other._impl);
        return *this; 
	};
	UnitVisible::~UnitVisible() = default;

	/********************
	Check setup
	********************/

	void UnitVisible::check_setup() const {
		_impl->check_setup();
	};

	/********************
	Location
	********************/

	int UnitVisible::dim() const {
		return _impl->dim();
	};
	int UnitVisible::x() const {
		return _impl->x();
	};
	int UnitVisible::y() const {
		return _impl->y();
	};
	int UnitVisible::z() const {
		return _impl->z();
	};
	bool UnitVisible::less_than(const UnitVisible &other) const {
		return _impl->less_than(other);
	};

	/********************
	Finish setup in lattice
	********************/

	// Add possible species
	void UnitVisible::add_possible_species(Sptr species) {
		_impl->add_possible_species(species);
	};
	void UnitVisible::set_possible_species(std::vector<Sptr> species) {
		_impl->set_possible_species(species);
	};
	const std::vector<Sptr>& UnitVisible::get_possible_species() const {
		return _impl->get_possible_species();
	};

	// Bias dict
	void UnitVisible::set_bias_dict(std::shared_ptr<BiasDict> bias_dict) {
		_impl->set_bias_dict(bias_dict);
	};
	const std::shared_ptr<BiasDict>& UnitVisible::get_bias_dict() const {
		return _impl->get_bias_dict();
	};

	// Add connection
	void UnitVisible::add_conn(ConnVV *conn, int idx_of_me) {
		_impl->add_conn(conn,idx_of_me);
	};
	std::vector<std::pair<ConnVV*,int>> UnitVisible::get_conns_vv() const {
		return _impl->get_conns_vv();
	};

	void UnitVisible::add_conn(ConnVVV *conn, int idx_of_me) {
		_impl->add_conn(conn,idx_of_me);
	};
	std::vector<std::pair<ConnVVV*,int>> UnitVisible::get_conns_vvv() const {
		return _impl->get_conns_vvv();
	};

	void UnitVisible::add_conn(ConnVH *conn) {
		_impl->add_conn(conn);
	};

	/********************
	Get probability
	********************/

	// Check mode
	bool UnitVisible::check_is_b_mode() const {
		return _impl->check_is_b_mode();
	};

	// Flip between the two modes
	void UnitVisible::set_b_mode(bool flag) {
		_impl->set_b_mode(flag);
	};

	// Binary
	Sptr UnitVisible::get_b_mode_species() const { // nullptr for empty
		return _impl->get_b_mode_species();
	};
	void UnitVisible::set_b_mode_species(Sptr sp) {
		_impl->set_b_mode_species(sp);
	};
	void UnitVisible::set_b_mode_species(std::string sp) {
		_impl->set_b_mode_species(sp);
	};
	void UnitVisible::set_b_mode_empty() {
		_impl->set_b_mode_empty();
	};
	bool UnitVisible::check_b_mode_is_empty() const {
		return _impl->check_b_mode_is_empty();
	};

	// Probabilistic
	double UnitVisible::get_p_mode_prob(Sptr sp) const { // nullptr for empty
		return _impl->get_p_mode_prob(sp);
	};
	const std::unordered_map<Sptr, double>& UnitVisible::get_p_mode_probs() const {
		return _impl->get_p_mode_probs();
	};
	void UnitVisible::set_p_mode_prob(Sptr sp, double prob) {
		_impl->set_p_mode_prob(sp,prob);
	};
	void UnitVisible::set_p_mode_prob(std::string sp, double prob) {
		_impl->set_p_mode_prob(sp,prob);
	};
	void UnitVisible::set_p_mode_empty() {
		_impl->set_p_mode_empty();
	};
	bool UnitVisible::check_p_mode_is_empty() const {
		return _impl->check_p_mode_is_empty();
	};

	// Convert between the two
	void UnitVisible::convert_b_to_p_mode() {
		_impl->convert_b_to_p_mode();
	};
	void UnitVisible::convert_p_to_b_mode() {
		_impl->convert_p_to_b_mode();
	};

	/********************
	Get moment
	********************/

	double UnitVisible::get_moment(std::string ixn_param_name, bool binary) const {
		return _impl->get_moment(ixn_param_name,binary);
	};

	/********************
	Get activation
	********************/

	double UnitVisible::get_act_for_species_at_timepoint(Sptr &species, int timepoint) const {
		return _impl->get_act_for_species_at_timepoint(species,timepoint);
	};

	/********************
	Sample
	********************/

	void UnitVisible::sample_at_timepoint(int timepoint, bool binary) {
		_impl->sample_at_timepoint(timepoint, binary);
	};
};