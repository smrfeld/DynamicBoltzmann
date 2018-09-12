#ifndef STRING_H
#define STRING_H
#include <string>
#endif

#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	Struct for specifying a dimension
	****************************************/

	// Type of dimension
	enum DimType { H, J, K, W, B };

	class Dim {

	private:

		class Impl;
		std::unique_ptr<Impl> _impl;

	public:

		/********************
		Constructor
		********************/

		Dim(std::string name, DimType type, std::vector<std::string> basis_func_dims, double min, double max, int n, double init);
		Dim(std::string name, DimType type, std::string species, std::vector<std::string> basis_func_dims, double min, double max, int n, double init);
		Dim(std::string name, DimType type, std::vector<std::string> species, std::vector<std::string> basis_func_dims, double min, double max, int n, double init);
		Dim(std::string name, DimType type, std::vector<std::vector<std::string>> species, std::vector<std::string> basis_func_dims, double min, double max, int n, double init);
		Dim(const Dim& other);
		Dim(Dim&& other);
		Dim& operator=(Dim other);
		~Dim();

		/********************
		Getters
		********************/

		// Name
		std::string name() const;

		// Type
		DimType type() const;

		// Does it apply to any species?
		bool any_species() const;

		// Basis func dims
		std::vector<std::string> basis_func_dims() const;

		// Min/max/n/init
		double min() const;
		double max() const;
		double n() const;
		double init() const;

		// Get species
		std::vector<std::string> get_species_h() const;
		std::vector<std::string> get_species_b() const;
		std::vector<std::vector<std::string>> get_species_J() const;
		std::vector<std::vector<std::string>> get_species_K() const;
		std::vector<std::vector<std::string>> get_species_W() const;

		/********************
		Setters
		********************/

		// Add basis func dimension
		void add_basis_func_dim(std::string dim);

		// Add species
		void add_species_h(std::string species);
		void add_species_b(std::string species);
		void add_species_J(std::string species1, std::string species2);
		void add_species_K(std::string species1, std::string species2, std::string species3);
		void add_species_W(std::string species_visible, std::string species_hidden);
	};

};