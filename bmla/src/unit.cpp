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

/************************************
* Namespace for bmla
************************************/

namespace bmla {

	/****************************************
	Class to hold a lattice site
	***************************************/

	class Unit::Impl {

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
		Impl(int x, int y);
		Impl(int x, int y, int z);
		Impl(int x, std::vector<Sptr> species_possible);
		Impl(int x, int y, std::vector<Sptr> species_possible);
		Impl(int x, int y, int z, std::vector<Sptr> species_possible);
		Impl(const Impl& other);
		Impl(Impl&& other);
		Impl& operator=(const Impl& other);
		Impl& operator=(Impl&& other);
		~Impl();

		/********************
		Verbose
		********************/

		std::string print_str() const;

		/********************
		Location
		********************/

		int dim() const;
		int x() const;
		int y() const;
		int z() const;
		bool less_than(const Unit &other) const;

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
		bool check_is_b_mode_species(Sptr sp) const;
		void set_b_mode_species(std::string sp);
		void set_b_mode_empty();
		bool check_b_mode_is_empty() const;
		void set_b_mode_random();

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
		Sample
		********************/

		void sample_given_activations(const std::vector<double>& activations, bool binary=true);
	};









































	/****************************************
	Class to hold a lattice Unit
	****************************************/

	/********************
	Constructor
	********************/

	Unit::Impl::Impl(int x) : Impl(x,0,0) { 
		_dim = 1;
	};
	Unit::Impl::Impl(int x, int y) : Impl(x,y,0) {
		_dim = 2;
	};
	Unit::Impl::Impl(int x, int y, int z) {
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
	Unit::Impl::Impl(int x, std::vector<Sptr> species_possible) : Impl(x,0,0,species_possible) { _dim=1; };
	Unit::Impl::Impl(int x, int y, std::vector<Sptr> species_possible) : Impl(x,y,0,species_possible) { _dim=2; };
	Unit::Impl::Impl(int x, int y, int z, std::vector<Sptr> species_possible) : Impl(x,y,z) {
		set_possible_species(species_possible);
	};	
	Unit::Impl::Impl(const Impl& other) {
		_copy(other);
	};
	Unit::Impl::Impl(Impl&& other) {
		_copy(other);
		other._reset();
	};
	Unit::Impl& Unit::Impl::operator=(const Impl& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	Unit::Impl& Unit::Impl::operator=(Impl&& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
	};
	Unit::Impl::~Impl() {
		_clean_up();
	};
	void Unit::Impl::_clean_up() {
		// Nothing....
	};
	void Unit::Impl::_reset() {
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

		_sampling_prob = 0.0;
		_sampling_rand = 0.0;
		_sampling_tot = 0.0;
		_sampling_i_chosen = 0;
		_sampling_props.clear();
		_sampling_probs.clear();
	};
	void Unit::Impl::_copy(const Impl& other) {
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

	void Unit::Impl::_sample_prop_vec() {
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

	std::string Unit::Impl::print_str() const {
		std::ostringstream s;
		if (_dim == 1) {
			s << _x << " ";
		} else if (_dim == 2) {
			s << _x << " " << _y << " ";
		} else if (_dim == 3) {
			s << _x << " " << _y << " " << _z << " ";
		};

		if (_b_mode_flag) {
			// Binary
			if (!_b_mode_sp) {
				s << "empty";
			} else {
				s << _b_mode_sp->get_name();
			};
		} else {
			// Prob
			s << "(empty," << _p_mode_prob_empty << ") ";
			for (auto &pr: _p_mode_probs_str) {
				s << "(" << pr.first << "," << *pr.second << ") ";
			};
		};	
		return s.str();
	};

	/********************
	Location
	********************/

	int Unit::Impl::dim() const {
		return _dim;
	};
	int Unit::Impl::x() const {
		return _x;
	};
	int Unit::Impl::y() const {
		return _y;
	};
	int Unit::Impl::z() const {
		return _z;
	};
	bool Unit::Impl::less_than(const Unit &other) const {
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
	void Unit::Impl::add_possible_species(Sptr species) {
		// Check that species is not already added
		auto it = std::find(_sp_possible.begin(), _sp_possible.end(), species);
		if (it != _sp_possible.end()) {
			std::cerr << ">>> Error: Unit::Impl::add_possible_species <<< species: " << species->get_name() << " already added." << std::endl;
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
	void Unit::Impl::set_possible_species(std::vector<Sptr> species) {
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
	const std::vector<Sptr>& Unit::Impl::get_possible_species() const {
		return _sp_possible;
	};

	// Bias dict
	void Unit::Impl::set_bias_dict(std::shared_ptr<BiasDict> bias_dict) {
		_bias_dict = bias_dict;
	};
	const std::shared_ptr<BiasDict>& Unit::Impl::get_bias_dict() const {
		return _bias_dict;
	};

	/********************
	Get probability
	********************/

	// Check mode
	bool Unit::Impl::check_is_b_mode() const {
		if (_b_mode_flag) {
			return true;
		} else {
			return false;
		};
	};

	// Flip between the two modes
	void Unit::Impl::set_b_mode(bool flag) {
		_b_mode_flag = flag;
	};

	// Binary
	Sptr Unit::Impl::get_b_mode_species() const {
		return _b_mode_sp;
	};
	void Unit::Impl::set_b_mode_species(Sptr sp) {
		_b_mode_sp = sp;
	};
	bool Unit::Impl::check_is_b_mode_species(Sptr sp) const {
		if (_b_mode_sp == sp) {
			return true;
		} else {
			return false;
		};
	};
	void Unit::Impl::set_b_mode_species(std::string sp) {
		_b_mode_sp = _sp_str_map[sp];
	};
	void Unit::Impl::set_b_mode_empty() {
		_b_mode_sp = nullptr;
	};
	bool Unit::Impl::check_b_mode_is_empty() const {
		if (_b_mode_sp) {
			return false;
		} else {
			return true;
		};
	};
	void Unit::Impl::set_b_mode_random() {
		int r = randI(0,_sp_possible.size());
		if (r == _sp_possible.size()) {
			set_b_mode_empty();
		} else {
			set_b_mode_species(_sp_possible[r]);
		};
	};

	// Probabilistic
	double Unit::Impl::get_p_mode_prob(Sptr sp) const {
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
	const std::unordered_map<Sptr, double>& Unit::Impl::get_p_mode_probs() const {
		return _p_mode_probs;
	};
	void Unit::Impl::set_p_mode_prob(Sptr sp, double prob) {
		// nullptr for prob of empty
		if (sp == nullptr) {
			_p_mode_prob_empty = prob;
			return;
		};

		_p_mode_probs[sp] = prob;
	};
	void Unit::Impl::set_p_mode_prob(std::string sp, double prob) {
		if (sp == "") {
			_p_mode_prob_empty = prob;
			return;
		};
		
		*(_p_mode_probs_str[sp]) = prob;
	};
	void Unit::Impl::set_p_mode_empty() {
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
	bool Unit::Impl::check_p_mode_is_empty() const {
		if (_p_mode_prob_empty == 1.0) {
			return true;
		} else {
			return false;
		};
	};

	// Convert between the two
	void Unit::Impl::convert_b_to_p_mode() {
		// Clear current
		set_p_mode_empty();

		// Occupancy
		if (_b_mode_sp) {
			set_p_mode_prob(_b_mode_sp,1.0);
			_p_mode_prob_empty = 0.0;
		};

		// Binary mode = false
		_b_mode_flag = false;
	};
	void Unit::Impl::convert_p_to_b_mode() {

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

	double Unit::Impl::get_moment(std::string ixn_param_name, bool binary) const {
		if (!_bias_dict) {
			std::cerr << ">>> Error: Unit::Impl::get_moment <<< no bias dict exists on this unit" << std::endl;
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
	Sample
	********************/

	void Unit::Impl::sample_given_activations(const std::vector<double>& activations, bool binary) {
		// 0 element of props is always 0
		// _sampling_props.push_back(0.0);

		// Empty = 1
		// 1st element of props = 0th element of probs is always 1
		// _sampling_props.push_back(1.0);
		// _sampling_probs.push_back(1.0);

		// Commit probs
		if (binary) {

			// Go through all possible species this could be, calculate propensities
			for (auto i=0; i<activations.size(); i++) {
				
				_sampling_prob = exp(activations[i]);

				// Append prop
				_sampling_props[i+2] = _sampling_props[i+1] + _sampling_prob;
			};

			/*
			std::cout << "Propensities: ";
			for (auto const &x: _sampling_props) {
				std::cout << x << " ";
			};
			std::cout << std::endl;
			*/

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
			_sampling_tot = 1.0;
			for (auto i=0; i<activations.size(); i++) {
				
				_sampling_prob = exp(activations[i]);

				// Append prob
				_sampling_probs[i+1] = _sampling_prob;
				_sampling_tot += _sampling_prob;
			};

			/*
			std::cout << "Sampling probs: " << std::flush;
			for (auto &x: _sampling_probs) {
				std::cout << x / _sampling_tot << " ";
			};
			std::cout << std::endl;
			*/
			
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

	Unit::Unit(int x) : _impl(new Impl(x)) {};
	Unit::Unit(int x, std::vector<Sptr> species_possible) : _impl(new Impl(x,species_possible)) {};
	Unit::Unit(int x, int y) : _impl(new Impl(x,y)) {};
	Unit::Unit(int x, int y, std::vector<Sptr> species_possible) : _impl(new Impl(x,y,species_possible)) {};
	Unit::Unit(int x, int y, int z) : _impl(new Impl(x,y,z)) {};
	Unit::Unit(int x, int y, int z, std::vector<Sptr> species_possible) : _impl(new Impl(x,y,z,species_possible)) {};
	Unit::Unit(const Unit& other) : _impl(new Impl(*other._impl)) {};
	Unit::Unit(Unit&& other) : _impl(std::move(other._impl)) {};
	Unit& Unit::operator=(const Unit& other) {
        _impl.reset( new Impl( *other._impl ) );
        return *this; 
	};
	Unit& Unit::operator=(Unit&& other) {
        _impl = std::move(other._impl);
        return *this; 
	};
	Unit::~Unit() = default;

	/********************
	Check setup
	********************/

	void Unit::print() const {
		std::cout << _impl->print_str() << std::endl;
	};
	std::string Unit::print_str() const {
		return _impl->print_str();
	};

	/********************
	Location
	********************/

	int Unit::dim() const {
		return _impl->dim();
	};
	int Unit::x() const {
		return _impl->x();
	};
	int Unit::y() const {
		return _impl->y();
	};
	int Unit::z() const {
		return _impl->z();
	};
	bool Unit::less_than(const Unit &other) const {
		return _impl->less_than(other);
	};

	/********************
	Finish setup in lattice
	********************/

	// Add possible species
	void Unit::add_possible_species(Sptr species) {
		_impl->add_possible_species(species);
	};
	void Unit::set_possible_species(std::vector<Sptr> species) {
		_impl->set_possible_species(species);
	};
	const std::vector<Sptr>& Unit::get_possible_species() const {
		return _impl->get_possible_species();
	};

	// Bias dict
	void Unit::set_bias_dict(std::shared_ptr<BiasDict> bias_dict) {
		_impl->set_bias_dict(bias_dict);
	};
	const std::shared_ptr<BiasDict>& Unit::get_bias_dict() const {
		return _impl->get_bias_dict();
	};

	/********************
	Get probability
	********************/

	// Check mode
	bool Unit::check_is_b_mode() const {
		return _impl->check_is_b_mode();
	};

	// Flip between the two modes
	void Unit::set_b_mode(bool flag) {
		_impl->set_b_mode(flag);
	};

	// Binary
	Sptr Unit::get_b_mode_species() const { // nullptr for empty
		return _impl->get_b_mode_species();
	};
	void Unit::set_b_mode_species(Sptr sp) {
		_impl->set_b_mode_species(sp);
	};
	bool Unit::check_is_b_mode_species(Sptr sp) const {
		return _impl->check_is_b_mode_species(sp);
	};
	void Unit::set_b_mode_species(std::string sp) {
		_impl->set_b_mode_species(sp);
	};
	void Unit::set_b_mode_empty() {
		_impl->set_b_mode_empty();
	};
	bool Unit::check_b_mode_is_empty() const {
		return _impl->check_b_mode_is_empty();
	};
	void Unit::set_b_mode_random() {
		return _impl->set_b_mode_random();
	};

	// Probabilistic
	double Unit::get_p_mode_prob(Sptr sp) const { // nullptr for empty
		return _impl->get_p_mode_prob(sp);
	};
	const std::unordered_map<Sptr, double>& Unit::get_p_mode_probs() const {
		return _impl->get_p_mode_probs();
	};
	void Unit::set_p_mode_prob(Sptr sp, double prob) {
		_impl->set_p_mode_prob(sp,prob);
	};
	void Unit::set_p_mode_prob(std::string sp, double prob) {
		_impl->set_p_mode_prob(sp,prob);
	};
	void Unit::set_p_mode_empty() {
		_impl->set_p_mode_empty();
	};
	bool Unit::check_p_mode_is_empty() const {
		return _impl->check_p_mode_is_empty();
	};

	// Convert between the two
	void Unit::convert_b_to_p_mode() {
		_impl->convert_b_to_p_mode();
	};
	void Unit::convert_p_to_b_mode() {
		_impl->convert_p_to_b_mode();
	};

	/********************
	Get moment
	********************/

	double Unit::get_moment(std::string ixn_param_name, bool binary) const {
		return _impl->get_moment(ixn_param_name,binary);
	};

	/********************
	Sample
	********************/

	void Unit::sample_given_activations(const std::vector<double>& activations, bool binary) {
		_impl->sample_given_activations(activations, binary);
	};












































	/****************************************
	Class to hold a lattice site
	***************************************/

	class UnitVisible::Impl {

	private:

		// Connections
		std::vector<std::pair<ConnVV*,int>> _conns_vv;
		std::vector<std::pair<ConnVVV*,int>> _conns_vvv;
		std::vector<ConnVH*> _conns_vh;

		// Activations
		std::vector<double> _activations;

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const Impl& other);

	public:

		/********************
		Constructor
		********************/

		Impl();
		Impl(const Impl& other);
		Impl(Impl&& other);
		Impl& operator=(const Impl& other);
		Impl& operator=(Impl&& other);
		~Impl();

		/********************
		Add connection
		********************/

		void add_conn(ConnVV *conn, int idx_of_me);
		const std::vector<std::pair<ConnVV*,int>>& get_conns_vv() const;

		void add_conn(ConnVVV *conn, int idx_of_me);
		const std::vector<std::pair<ConnVVV*,int>>& get_conns_vvv() const;
		
		void add_conn(ConnVH *conn);
		const std::vector<ConnVH*>& get_conns_vh() const;

		/********************
		Get activation
		********************/

		double get_activation_for_species(const std::shared_ptr<BiasDict>& bias_dict, const Sptr &species) const;

		/********************
		Sample
		********************/

		void form_activations_vector(const std::shared_ptr<BiasDict>& bias_dict, const std::vector<Sptr>& species);
		const std::vector<double>& get_activations_vector() const;
	};



 








































	/****************************************
	Class to hold a lattice UnitVisible
	****************************************/

	/********************
	Constructor
	********************/

	UnitVisible::Impl::Impl() {};
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
		_conns_vv.clear();
		_conns_vvv.clear();
		_conns_vh.clear();
		_activations.clear();
	};
	void UnitVisible::Impl::_copy(const Impl& other) {
		_conns_vv = other._conns_vv;
		_conns_vvv = other._conns_vvv;
		_conns_vh = other._conns_vh;
		_activations = other._activations;
	};

	/********************
	Add connection
	********************/

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
	const std::vector<std::pair<ConnVV*,int>>& UnitVisible::Impl::get_conns_vv() const {
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
	const std::vector<std::pair<ConnVVV*,int>>& UnitVisible::Impl::get_conns_vvv() const {
		return _conns_vvv;
	};

	void UnitVisible::Impl::add_conn(ConnVH *conn) {
		_conns_vh.push_back(conn);
	};
	const std::vector<ConnVH*>& UnitVisible::Impl::get_conns_vh() const {
		return _conns_vh;
	};

	/********************
	Get activation
	********************/

	double UnitVisible::Impl::get_activation_for_species(const std::shared_ptr<BiasDict>& bias_dict, const Sptr &species) const {
		double act = 0.0;

		// Bias
		if (bias_dict) {
			act += bias_dict->get_ixn(species);
		};

		// VV conns
		for (auto const &conn_vv: _conns_vv) {
			act += conn_vv.first->get_act_for_species_at_unit(species,conn_vv.second);
		};

		// VVV conns
		for (auto const &conn_vvv: _conns_vvv) {
			act += conn_vvv.first->get_act_for_species_at_unit(species,conn_vvv.second);
		};

		// VH conns
		for (auto const &conn_vh: _conns_vh) {
			act += conn_vh->get_act_for_species_at_unit_v(species);
		};

		return act;
	};

	/********************
	Sample
	********************/

	void UnitVisible::Impl::form_activations_vector(const std::shared_ptr<BiasDict>& bias_dict, const std::vector<Sptr>& species) {
		// Check size
		if (_activations.size() != species.size()) {
			_activations.clear();
			for (auto i=0; i<species.size(); i++) {
				_activations.push_back(0.0);
			};
		};

		// Form the vector of activations
		int i=0;
		for (auto const &sp: species) {
			_activations[i] = get_activation_for_species(bias_dict,sp);
			i++;
		};
	};

	const std::vector<double>& UnitVisible::Impl::get_activations_vector() const {
		return _activations;
	};












































	/****************************************
	Forwards
	****************************************/

	/********************
	Constructor
	********************/

	UnitVisible::UnitVisible(int x) : Unit(x), _impl(new Impl()) {};
	UnitVisible::UnitVisible(int x, int y) : Unit(x,y), _impl(new Impl()) {};
	UnitVisible::UnitVisible(int x, int y, int z) : Unit(x,y,z), _impl(new Impl()) {};

	UnitVisible::UnitVisible(int x, std::vector<Sptr> species_possible) : Unit(x,species_possible), _impl(new Impl()) {};
	UnitVisible::UnitVisible(int x, int y, std::vector<Sptr> species_possible) : Unit(x,y,species_possible), _impl(new Impl()) {};
	UnitVisible::UnitVisible(int x, int y, int z, std::vector<Sptr> species_possible) : Unit(x,y,z,species_possible), _impl(new Impl()) {};

	UnitVisible::UnitVisible(const UnitVisible& other) : Unit(other), _impl(new Impl(*other._impl)) {};
	UnitVisible::UnitVisible(UnitVisible&& other) : Unit(other), _impl(std::move(other._impl)) {};
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
	Verbose
	********************/

	void UnitVisible::print() const {
		std::cout << print_str() << std::endl;
	};
	std::string UnitVisible::print_str() const {
		return "[Visible] " + Unit::print_str();
	};

	/********************
	Add connection
	********************/

	void UnitVisible::add_conn(ConnVV *conn, int idx_of_me) {
		_impl->add_conn(conn,idx_of_me);
	};
	const std::vector<std::pair<ConnVV*,int>>& UnitVisible::get_conns_vv() const {
		return _impl->get_conns_vv();
	};

	void UnitVisible::add_conn(ConnVVV *conn, int idx_of_me) {
		_impl->add_conn(conn,idx_of_me);
	};
	const std::vector<std::pair<ConnVVV*,int>>& UnitVisible::get_conns_vvv() const {
		return _impl->get_conns_vvv();
	};

	void UnitVisible::add_conn(ConnVH *conn) {
		_impl->add_conn(conn);
	};
	const std::vector<ConnVH*>& UnitVisible::get_conns_vh() const {
		return _impl->get_conns_vh();
	};

	/********************
	Get activation
	********************/

	double UnitVisible::get_activation_for_species(Sptr &species) const {
		return _impl->get_activation_for_species(get_bias_dict(),species);
	};

	/********************
	Sample
	********************/

	void UnitVisible::sample(bool binary) {
		_impl->form_activations_vector(get_bias_dict(),get_possible_species());

		// Parent sampler
		sample_given_activations(_impl->get_activations_vector(),binary);
	};
















































	/****************************************
	Class to hold a lattice site
	***************************************/

	class UnitHidden::Impl {

	private:

		// Layer
		int _layer;

		// Connections
		std::vector<ConnVH*> _conns_vh;
		// Layer -> conns -> (conn,idx of me in conn)
		std::map<int,std::vector<std::pair<ConnHH*,int>>> _conns_hh;

		// Activations
		std::vector<double> _activations;

		// Constructor helpers
		void _clean_up();
		void _reset();
		void _copy(const Impl& other);

	public:

		/********************
		Constructor
		********************/

		Impl(int layer);
		Impl(const Impl& other);
		Impl(Impl&& other);
		Impl& operator=(const Impl& other);
		Impl& operator=(Impl&& other);
		~Impl();

		/********************
		Location
		********************/

		int layer() const;

		/********************
		Add connection
		********************/

		// Visible-hidden conns
		void add_conn(ConnVH *conn);
		const std::vector<ConnVH*>& get_conns_vh() const;

		// Hidden-hidden conns
		void add_conn(ConnHH *conn, int idx_of_me, int from_layer);
		const std::map<int,std::vector<std::pair<ConnHH*,int>>>& get_conns_hh() const;

		/********************
		Get activation
		********************/

		double get_activation_for_species(const std::shared_ptr<BiasDict>& bias_dict, const Sptr &species, int given_layer) const;

		/********************
		Sample
		********************/

		void form_activations_vector(const std::shared_ptr<BiasDict>& bias_dict, const std::vector<Sptr>& species, int given_layer);
		const std::vector<double>& get_activations_vector() const;
	};



 








































	/****************************************
	Class to hold a lattice UnitHidden
	****************************************/

	/********************
	Constructor
	********************/

	UnitHidden::Impl::Impl(int layer) {
		_layer = layer;
	};
	UnitHidden::Impl::Impl(const Impl& other) {
		_copy(other);
	};
	UnitHidden::Impl::Impl(Impl&& other) {
		_copy(other);
		other._reset();
	};
	UnitHidden::Impl& UnitHidden::Impl::operator=(const Impl& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	UnitHidden::Impl& UnitHidden::Impl::operator=(Impl&& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
	};
	UnitHidden::Impl::~Impl() {
		_clean_up();
	};
	void UnitHidden::Impl::_clean_up() {
		// Nothing....
	};
	void UnitHidden::Impl::_reset() {
		_layer = 0;
		_activations.clear();
		_conns_vh.clear();
		_conns_hh.clear();
	};
	void UnitHidden::Impl::_copy(const Impl& other) {
		_layer = other._layer;
		_activations = other._activations;
		_conns_vh = other._conns_vh;
		_conns_hh = other._conns_hh;
	};

	/********************
	Location
	********************/

	int UnitHidden::Impl::layer() const {
		return _layer;
	};

	/********************
	Add connection
	********************/

	void UnitHidden::Impl::add_conn(ConnVH *conn) {
		_conns_vh.push_back(conn);
	};
	const std::vector<ConnVH*>& UnitHidden::Impl::get_conns_vh() const {
		return _conns_vh;
	};

	// Hidden-hidden conns
	void UnitHidden::Impl::add_conn(ConnHH *conn, int idx_of_me, int from_layer) {

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
				std::cerr << ">>> Error: UnitHidden::Impl::add_conn <<< connection already exists!" << std::endl;
				exit(EXIT_FAILURE);	
			} else {
				// Add
				_conns_hh[from_layer].push_back(pr);
			};
		};
	};
	const std::map<int,std::vector<std::pair<ConnHH*,int>>>& UnitHidden::Impl::get_conns_hh() const {
		return _conns_hh;
	};

	/********************
	Get activation
	********************/

	double UnitHidden::Impl::get_activation_for_species(const std::shared_ptr<BiasDict>& bias_dict, const Sptr &species, int given_layer) const {
		double act = 0.0;

		// Bias
		if (bias_dict) {
			act += bias_dict->get_ixn(species);
		};

		// VH conns
		if (given_layer == 0) {
			for (auto const &conn_vh: _conns_vh) {
				act += conn_vh->get_act_for_species_at_unit_h(species);
			};
		};

		// HH conns
		auto layer = _conns_hh.find(given_layer);
		if (layer != _conns_hh.end()) {
			// Go through conns
			for (auto const &conn_hh: layer->second) {
				act += conn_hh.first->get_act_for_species_at_unit(species,conn_hh.second);
			};
		};

		return act;
	};

	/********************
	Sample
	********************/

	void UnitHidden::Impl::form_activations_vector(const std::shared_ptr<BiasDict>& bias_dict, const std::vector<Sptr>& species, int given_layer) {
		// Check size
		if (_activations.size() != species.size()) {
			_activations.clear();
			for (auto i=0; i<species.size(); i++) {
				_activations.push_back(0.0);
			};
		};

		// Form the vector of activations
		for (auto i=0; i<species.size(); i++) {
			_activations[i] = get_activation_for_species(bias_dict,species[i], given_layer);
		};
	};

	const std::vector<double>& UnitHidden::Impl::get_activations_vector() const {
		return _activations;
	};





































	/****************************************
	Forwards
	****************************************/

	/********************
	Constructor
	********************/

	UnitHidden::UnitHidden(int layer, int idx) : Unit(idx), _impl(new Impl(layer)) {};

	UnitHidden::UnitHidden(int layer, int idx, std::vector<Sptr> species_possible) : Unit(idx,species_possible), _impl(new Impl(layer)) {};

	UnitHidden::UnitHidden(const UnitHidden& other) : Unit(other), _impl(new Impl(*other._impl)) {};
	UnitHidden::UnitHidden(UnitHidden&& other) : Unit(other), _impl(std::move(other._impl)) {};
	UnitHidden& UnitHidden::operator=(const UnitHidden& other) {
        _impl.reset( new Impl( *other._impl ) );
        return *this; 
	};
	UnitHidden& UnitHidden::operator=(UnitHidden&& other) {
        _impl = std::move(other._impl);
        return *this; 
	};
	UnitHidden::~UnitHidden() = default;
	
	/********************
	Verbose
	********************/

	void UnitHidden::print() const {
		std::cout << print_str() << std::endl;
	};
	std::string UnitHidden::print_str() const {
		std::ostringstream s;
		s << "[Hidden: layer: " << layer() << "] " << Unit::print_str();
		return s.str();
	};

	/********************
	Location
	********************/

	int UnitHidden::layer() const {
		return _impl->layer();
	};
	int UnitHidden::idx() const {
		return Unit::x();
	};

	/********************
	Add connection
	********************/

	void UnitHidden::add_conn(ConnVH *conn) {
		_impl->add_conn(conn);
	};
	const std::vector<ConnVH*>& UnitHidden::get_conns_vh() const {
		return _impl->get_conns_vh();
	};

	// Hidden-hidden conns
	void UnitHidden::add_conn(ConnHH *conn, int idx_of_me, int from_layer) {
		_impl->add_conn(conn,idx_of_me,from_layer);
	};
	const std::map<int,std::vector<std::pair<ConnHH*,int>>>& UnitHidden::get_conns_hh() const {
		return _impl->get_conns_hh();
	};

	/********************
	Get activation
	********************/

	double UnitHidden::get_activation_for_species(Sptr &species, int given_layer) const {
		return _impl->get_activation_for_species(get_bias_dict(),species, given_layer);
	};

	/********************
	Sample
	********************/

	void UnitHidden::sample(int given_layer, bool binary) {
		_impl->form_activations_vector(get_bias_dict(),get_possible_species(), given_layer);

		// Parent sampler
		sample_given_activations(_impl->get_activations_vector(),binary);
	};
};