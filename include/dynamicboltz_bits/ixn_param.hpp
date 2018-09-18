#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

#ifndef STRING_H
#define STRING_H
#include <string>
#endif

#ifndef FWDS_SPECIES_H
#define FWDS_SPECIES_H
#include "fwds/fwds_species.hpp"
#endif

#ifndef FWDS_IXN_PARAM_H
#define FWDS_IXN_PARAM_H
#include "fwds/fwds_ixn_param.hpp"
#endif

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	// Forwards
	class DiffEqRHS;
	class Moment;

	/****************************************
	Interaction parameter
	****************************************/

	// Enumeration of type of dimension
	enum class IxnParamType: unsigned int { H, J, K, W, B };

	class IxnParam {

	private:

		class Impl;
		std::unique_ptr<Impl> _impl;

	public:

		/********************
		Constructor
		********************/

		IxnParam(std::string name, IxnParamType type, double init_cond); 
		IxnParam(const IxnParam& other);
		IxnParam(IxnParam&& other);
		IxnParam& operator=(const IxnParam& other);
		IxnParam& operator=(IxnParam&& other);
		~IxnParam();

		/********************
		Check setup
		********************/

		void check_setup() const;

		/********************
		Add species
		********************/

		void add_species_h(Sptr species);
		void add_species_b(Sptr species);
		void add_species_j(Sptr species_of_site_1, Sptr species_of_site_2);
		void add_species_k(Sptr species_of_site_1, Sptr species_of_site_2, Sptr species_of_site_3);
		void add_species_w(Sptr species_of_visible, Sptr species_of_hidden);

		/********************
		Get species
		********************/

		const std::vector<Sptr>& get_species_h() const;
		const std::vector<Sptr>& get_species_b() const;
		const std::vector<Sptr2>& get_species_j() const;
		const std::vector<Sptr3>& get_species_k() const;
		const std::vector<Sptr2>& get_species_w() const;

		/********************
		Timesteps
		********************/

		int get_no_timesteps() const;
		void set_no_timesteps(int no_timesteps);

		/********************
		Init cond
		********************/

		double get_init_cond() const;
		void set_init_cond(double init_cond);

		// void set_fixed_awake_moment(std::vector<double> vals);

		/********************
		Name, type
		********************/

		std::string get_name() const;

		IxnParamType get_type() const;

		/********************
		Value
		********************/

		double get_val_at_timepoint(int timepoint) const;

		/********************
		Diff eq
		********************/

		void set_diff_eq(std::shared_ptr<DiffEqRHS> diff_eq);

		void solve_diff_eq_at_timepoint_to_plus_one(int timepoint, double dt);

		/********************
		Moment
		********************/

		std::shared_ptr<Moment> get_moment() const;

		/********************
		Write into an ofstream
		********************/

		/*
		void write_vals(std::string dir, int idx_opt_step, std::vector<int> idxs, int n_t_traj) const;
		void write_moments(std::string dir, int idx_opt_step, std::vector<int> idxs, int n_t_traj) const;
		*/

	};

};

