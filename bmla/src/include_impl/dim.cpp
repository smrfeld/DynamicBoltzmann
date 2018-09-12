#include "../../include/bmla_bits/dim.hpp"

#include <iostream>

/************************************
* Namespace for bmla
************************************/

namespace bmla {

	/****************************************
	Dim - IMPLEMENTATION
	****************************************/

	class Dim::Impl {
	private:

		// Name
		std::string _name;

		// Type
		DimType _type;

		// Name of associated species
		bool _any_species;
		std::vector<std::string> _species;
		std::vector<std::vector<std::string> > _species_multiple;

		// Guess
		double _guess;

		// Constructor helpers
		void _clean_up();
		void _copy(const Impl& other);
		void _reset();
		void _shared_constructor(std::string name, DimType type, double guess);

	public:

		// Constructor
		Impl(std::string name, DimType type, double guess);
		Impl(std::string name, DimType type, std::string species, double guess);
		Impl(std::string name, DimType type, std::vector<std::string> species, double guess);
		Impl(std::string name, DimType type, std::vector<std::vector<std::string>> species, double guess);
		Impl(const Impl& other);
		Impl(Impl&& other);
	    Impl& operator=(const Impl& other);
	    Impl& operator=(Impl&& other);
		~Impl();

		/********************
		Getters
		********************/

		// Name
		std::string name() const;

		// Type
		DimType type() const;

		// Does it apply to all species?
		bool any_species() const;

		// Guess
		double guess() const;

		// Get species
		std::vector<std::string> get_species_h() const;
		std::vector<std::string> get_species_b() const;
		std::vector<std::vector<std::string>> get_species_J() const;
		std::vector<std::vector<std::string>> get_species_K() const;
		std::vector<std::vector<std::string>> get_species_W() const;
		std::vector<std::vector<std::string>> get_species_Q() const;

		/********************
		Setters
		********************/

		// Add species
		void add_species_h(std::string species);
		void add_species_b(std::string species);
		void add_species_J(std::string species1, std::string species2);
		void add_species_K(std::string species1, std::string species2, std::string species3);
		void add_species_W(std::string species_visible, std::string species_hidden);
		void add_species_Q(std::string species1, std::string species2, std::string species3, std::string species4);
	};

	/********************
	Constructor
	********************/

	Dim::Impl::Impl(std::string name, DimType type, double guess) {
		_shared_constructor(name,type,guess);
		// Yes any species
		_any_species = true;
	};
	Dim::Impl::Impl(std::string name, DimType type, std::string species, double guess) {
		_shared_constructor(name,type,guess);
		// Not any species
		_any_species = false;
		// Add
		if (_type == H) {
			add_species_h(species);
		} else if (_type == B) {
			add_species_b(species);
		};
	};
	Dim::Impl::Impl(std::string name, DimType type, std::vector<std::string> species, double guess) {
		_shared_constructor(name,type,guess);
		// Not any species
		_any_species = false;
		// Add
		if (_type == H) {
			for (auto s: species) {
				add_species_h(s);
			};
		} else if (_type == B) {
			for (auto s: species) {
				add_species_b(s);
			};
		} else if (type == J) {
			if (species.size() != 2) {
				std::cerr << "Error! must be 2 species for J" << std::endl;
				exit(EXIT_FAILURE);
			};
			add_species_J(species[0],species[1]);
		} else if (type == K) {
			if (species.size() != 3) {
				std::cerr << "Error! must be 3 species for K" << std::endl;
				exit(EXIT_FAILURE);
			};
			add_species_K(species[0],species[1],species[2]);
		} else if (type == W) {
			if (species.size() != 2) {
				std::cerr << "Error! must be 2 species for W" << std::endl;
				exit(EXIT_FAILURE);
			};
			add_species_W(species[0],species[1]);
 		} else if (type == Q) {
			if (species.size() != 4) {
				std::cerr << "Error! must be 4 species for Q" << std::endl;
				exit(EXIT_FAILURE);
			};
			add_species_Q(species[0],species[1],species[2],species[3]);
		};
	};
	Dim::Impl::Impl(std::string name, DimType type, std::vector<std::vector<std::string>> species, double guess) {
		// Check type
		if (type == B || type != H) {
			std::cerr << "Error! Only for J or K or W or Q." << std::endl;
			exit(EXIT_FAILURE);
		};
		_shared_constructor(name,type,guess);
		// Not any species
		_any_species = false;
		// Add
		if (type == J) {
			for (auto s_pair: species) {
				if (s_pair.size() != 2) {
					std::cerr << "Error! must be 2 species for J" << std::endl;
					exit(EXIT_FAILURE);
				};
				add_species_J(s_pair[0],s_pair[1]);
			};
		} else if (type == K) {
			for (auto s_triplet: species) {
				if (s_triplet.size() != 3) {
					std::cerr << "Error! must be 3 species for K" << std::endl;
					exit(EXIT_FAILURE);
				};
				add_species_K(s_triplet[0],s_triplet[1],s_triplet[2]);
			};
		} else if (type == W) {
			for (auto s_pair: species) {
				if (s_pair.size() != 2) {
					std::cerr << "Error! must be 2 species for W" << std::endl;
					exit(EXIT_FAILURE);
				};
				add_species_W(s_pair[0],s_pair[1]);
			};
		} else if (type == Q) {
			for (auto s_quartic: species) {
				if (s_quartic.size() != 4) {
					std::cerr << "Error! must be 4 species for Q" << std::endl;
					exit(EXIT_FAILURE);
				};
				add_species_Q(s_quartic[0],s_quartic[1],s_quartic[2],s_quartic[3]);
			};
		};
	};
	Dim::Impl::Impl(const Impl& other) {
		_copy(other);
	};
	Dim::Impl::Impl(Impl&& other) {
		_copy(other);
		other._reset();
	};
    Dim::Impl& Dim::Impl::operator=(const Impl& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    Dim::Impl& Dim::Impl::operator=(Impl&& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
    };
	Dim::Impl::~Impl()
	{
		_clean_up();
	};
	void Dim::Impl::_clean_up() {
		// Nothing...
	};
	void Dim::Impl::_copy(const Impl& other) {
		_name = other._name;
		_type = other._type;
		_any_species = other._any_species;
		_species = other._species;
		_species_multiple = other._species_multiple;
		_guess = other._guess;
	};
	void Dim::Impl::_reset() {
		_name = "";
		_any_species = false;
		_species.clear();
		_species_multiple.clear();
		_guess = 0.;
	};
	void Dim::Impl::_shared_constructor(std::string name, DimType type, double guess) {
		_name = name;
		_type = type;
		_guess = guess;
	};

	/********************
	Getters
	********************/

	// Name
	std::string Dim::Impl::name() const {
		return _name;
	};

	// Type
	DimType Dim::Impl::type() const {
		return _type;
	};

	// Does it apply to all species?
	bool Dim::Impl::any_species() const {
		return _any_species;
	};

	// Guess
	double Dim::Impl::guess() const {
		return _guess;
	};

	// Get species
	std::vector<std::string> Dim::Impl::get_species_h() const {
		if (_type != H) {
			std::cerr << "Error! Requested species but not of type H." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _species;
	};
	std::vector<std::string> Dim::Impl::get_species_b() const {
		if (_type != B) {
			std::cerr << "Error! Requested species but not of type B." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _species;
	};
	std::vector<std::vector<std::string>> Dim::Impl::get_species_J() const {
		if (_type != J) {
			std::cerr << "Error! Requested species but not of type J." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _species_multiple;
	};
	std::vector<std::vector<std::string>> Dim::Impl::get_species_K() const {
		if (_type != K) {
			std::cerr << "Error! Requested species but not of type K." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _species_multiple;
	};
	std::vector<std::vector<std::string>> Dim::Impl::get_species_W() const {
		if (_type != W) {
			std::cerr << "Error! Requested species but not of type W." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _species_multiple;
	};
	std::vector<std::vector<std::string>> Dim::Impl::get_species_Q() const {
		if (_type != Q) {
			std::cerr << "Error! Requested species but not of type Q." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _species_multiple;
	};

	/********************
	Setters
	********************/

	// Add species
	void Dim::Impl::add_species_h(std::string species) {
		// Check type
		if (_type != H) {
			std::cerr << "Error! Not H." << std::endl;
			exit(EXIT_FAILURE);
		};
		_any_species = false;
		_species.push_back(species);
	};
	void Dim::Impl::add_species_b(std::string species) {
		// Check type
		if (_type != B) {
			std::cerr << "Error! Not B." << std::endl;
			exit(EXIT_FAILURE);
		};
		_any_species = false;
		_species.push_back(species);
	};
	void Dim::Impl::add_species_J(std::string species1, std::string species2) {
		// Check type
		if (_type != J) {
			std::cerr << "Error! Not J." << std::endl;
			exit(EXIT_FAILURE);
		};
		_any_species = false;
		_species_multiple.push_back(std::vector<std::string>({species1,species2}));
	};
	void Dim::Impl::add_species_K(std::string species1, std::string species2, std::string species3) {
		// Check type
		if (_type != K) {
			std::cerr << "Error! Not K." << std::endl;
			exit(EXIT_FAILURE);
		};
		_any_species = false;
		_species_multiple.push_back(std::vector<std::string>({species1,species2,species3}));
	};
	void Dim::Impl::add_species_W(std::string species_visible, std::string species_hidden) {
		// Check type
		if (_type != W) {
			std::cerr << "Error! Not W." << std::endl;
			exit(EXIT_FAILURE);
		};
		_any_species = false;
		_species_multiple.push_back(std::vector<std::string>({species_visible,species_hidden}));
	};
	void Dim::Impl::add_species_Q(std::string species1, std::string species2, std::string species3, std::string species4) {
		// Check type
		if (_type != Q) {
			std::cerr << "Error! Not Q." << std::endl;
			exit(EXIT_FAILURE);
		};
		_any_species = false;
		_species_multiple.push_back(std::vector<std::string>({species1,species2,species3,species4}));
	};





















































	/****************************************
	Dim - Impl forwards
	****************************************/

	Dim::Dim(std::string name, DimType type, double guess) : _impl(new Impl(name,type,guess)) {};
	Dim::Dim(std::string name, DimType type, std::string species, double guess) : _impl(new Impl(name,type,species,guess)) {};
	Dim::Dim(std::string name, DimType type, std::vector<std::string> species, double guess) : _impl(new Impl(name,type,species,guess)) {};
	Dim::Dim(std::string name, DimType type, std::vector<std::vector<std::string>> species, double guess) : _impl(new Impl(name,type,species,guess)) {};
	Dim::Dim(const Dim& other) : _impl(new Impl(*other._impl)) {};
	Dim::Dim(Dim&& other) = default;
	Dim& Dim::operator=(Dim other) {
        _impl = std::move(other._impl);
        return *this; 
	};
	Dim::~Dim() = default;

	// Name
	std::string Dim::name() const {
		return _impl->name();
	};

	// Type
	DimType Dim::type() const {
		return _impl->type();
	};

	// Does it apply to all species?
	bool Dim::any_species() const {
		return _impl->any_species();
	};

	// Guess
	double Dim::guess() const {
		return _impl->guess();
	};

	// Get species
	std::vector<std::string> Dim::get_species_h() const {
		return _impl->get_species_h();
	};
	std::vector<std::string> Dim::get_species_b() const {
		return _impl->get_species_b();
	};
	std::vector<std::vector<std::string>> Dim::get_species_J() const {
		return _impl->get_species_J();
	};
	std::vector<std::vector<std::string>> Dim::get_species_K() const {
		return _impl->get_species_K();
	};
	std::vector<std::vector<std::string>> Dim::get_species_W() const {
		return _impl->get_species_W();
	};
	std::vector<std::vector<std::string>> Dim::get_species_Q() const {
		return _impl->get_species_Q();
	};

	/********************
	Setters
	********************/

	// Add species
	void Dim::add_species_h(std::string species) {
		_impl->add_species_h(species);
	};
	void Dim::add_species_b(std::string species) {
		_impl->add_species_b(species);
	};
	void Dim::add_species_J(std::string species1, std::string species2) {
		_impl->add_species_J(species1,species2);
	};
	void Dim::add_species_K(std::string species1, std::string species2, std::string species3) {
		_impl->add_species_K(species1,species2,species3);
	};
	void Dim::add_species_W(std::string species_visible, std::string species_hidden) {
		_impl->add_species_W(species_visible,species_hidden);
	};
	void Dim::add_species_Q(std::string species1, std::string species2, std::string species3, std::string species4) {
		_impl->add_species_Q(species1,species2,species3,species4);
	};
};


