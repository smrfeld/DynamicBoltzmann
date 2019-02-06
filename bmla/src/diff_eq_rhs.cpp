#include "../include/dblz_bits/diff_eq_rhs.hpp"

// Other headers
#include "../include/dblz_bits/general.hpp"
#include "../include/dblz_bits/ixn_param_traj.hpp"
#include "../include/dblz_bits/ixn_param.hpp"
#include "../include/dblz_bits/adjoint.hpp"

#include <iostream>
#include "math.h"
#include <fstream>
#include <iomanip>

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	Domain1D
	****************************************/

	Domain1D::Domain1D(ITptr ixn_param_traj, double delta, double zero) : Dimension1D(delta,zero) {
		_ixn_param_traj = ixn_param_traj;
	};
	Domain1D::Domain1D(const Domain1D& other) : Dimension1D(other) {
		_copy(other);
	};
	Domain1D::Domain1D(Domain1D&& other) : Dimension1D(std::move(other)) {
		_copy(other);
		other._reset();
	};
	Domain1D& Domain1D::operator=(const Domain1D& other) {
		if (this != &other)
		{
	        Dimension1D::operator=(other);
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	Domain1D& Domain1D::operator=(Domain1D&& other) {
		if (this != &other)
		{
	        Dimension1D::operator=(std::move(other));
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
	};
	Domain1D::~Domain1D() {
		_clean_up();
	};
	void Domain1D::_copy(const Domain1D& other)
	{
		_ixn_param_traj = other._ixn_param_traj;
	};
	void Domain1D::_reset()
	{
		_ixn_param_traj = nullptr;
	};
	void Domain1D::_clean_up() {
	};

	/********************
	Getters
	********************/

	std::string Domain1D::get_name() const {
		return _ixn_param_traj->get_name();
	};
	ITptr Domain1D::get_ixn_param_traj() const {
		return _ixn_param_traj;
	};












































	/****************************************
	DiffEqRHS
	****************************************/

	/********************
	Constructor
	********************/

	DiffEqRHS::DiffEqRHS(std::string name, ITptr parent_ixn_param_traj, std::vector<Domain1D*> domain) : q3c1::Grid(std::vector<q3c1::Dimension1D*>(domain.begin(),domain.end())) {
		_name = name;
		_domain = domain;
		_no_dims = _domain.size();
		_parent_ixn_param_traj = parent_ixn_param_traj;

		if (_no_dims == 0 || _no_dims > 3) {
			std::cerr << ">>> Error: DiffEqRHS::DiffEqRHS <<< Only dims 1,2,3 are supported" << std::endl;
			exit(EXIT_FAILURE);
		};
	};
	DiffEqRHS::DiffEqRHS(const DiffEqRHS& other) : Grid(other) {
		_copy(other);
	};
	DiffEqRHS::DiffEqRHS(DiffEqRHS&& other) : Grid(std::move(other)) {
		_move(other);
	};

	DiffEqRHS& DiffEqRHS::operator=(const DiffEqRHS& other) {
		if (this != &other)
		{
			_clean_up();

			Grid::operator=(other);
			_copy(other);
		};
		return *this;
	};
	DiffEqRHS& DiffEqRHS::operator=(DiffEqRHS&& other) {
		if (this != &other)
		{
			_clean_up();
			Grid::operator=(other);
			_move(other);
		};
		return *this;
	};

	DiffEqRHS::~DiffEqRHS() {
		_clean_up();
	};
	void DiffEqRHS::_copy(const DiffEqRHS& other) {
		_name = other._name;
		_no_dims = other._no_dims;
		_domain = other._domain;
		_parent_ixn_param_traj = other._parent_ixn_param_traj;
		_updates = other._updates;

        if (other._nesterov_y_s) {
			_nesterov_y_s = new std::map<q3c1::Vertex*,std::vector<double>>(*other._nesterov_y_s);
		};
		if (other._nesterov_y_sp1) {
			_nesterov_y_sp1 = new std::map<q3c1::Vertex*,std::vector<double>>(*other._nesterov_y_sp1);
		};
        if (other._adam_m) {
            _adam_m = new std::map<q3c1::Vertex*,std::vector<double>>(*other._adam_m);
        };
        if (other._adam_v) {
            _adam_v = new std::map<q3c1::Vertex*,std::vector<double>>(*other._adam_v);
        };
	};
	void DiffEqRHS::_move(DiffEqRHS& other) {
		_name = other._name;
		_no_dims = other._no_dims;
		_domain = other._domain;
		_parent_ixn_param_traj = other._parent_ixn_param_traj;
		_updates = other._updates;
        
		_nesterov_y_s = other._nesterov_y_s;
		_nesterov_y_sp1 = other._nesterov_y_sp1;

        _adam_m = other._adam_m;
        _adam_v = other._adam_v;
        
		// Reset other
		other._name = "";
		other._no_dims = 0;
		other._domain.clear();
		other._parent_ixn_param_traj = nullptr;
		other._updates.clear();

        other._nesterov_y_s = nullptr;
		other._nesterov_y_sp1 = nullptr;
        
        other._adam_m = nullptr;
        other._adam_v = nullptr;
	};

	void DiffEqRHS::_clean_up() {
		// Note: Domain is not owned => do not delete
		if (_nesterov_y_s) {
			delete _nesterov_y_s;
			_nesterov_y_s = nullptr;
		};
		if (_nesterov_y_sp1) {
			delete _nesterov_y_sp1;
			_nesterov_y_sp1 = nullptr;
		};
        if (_adam_m) {
            delete _adam_m;
            _adam_m = nullptr;
        };
        if (_adam_v) {
            delete _adam_v;
            _adam_v = nullptr;
        };
	};

	/********************
	Validate
	********************/

	// Validate setup
	void DiffEqRHS::validate_setup() const {
		std::cout << "--- Validate Basis func " << _name << " ---" << std::endl;
	};

	/********************
	Getters
	********************/

	std::string DiffEqRHS::get_name() const {
		return _name;
	};

	ITptr DiffEqRHS::get_parent_ixn_param_traj() const {
		return _parent_ixn_param_traj;
	};

	const std::vector<Domain1D*>& DiffEqRHS::get_domain() const {
		return _domain;
	};

	std::vector<double> DiffEqRHS::_form_abscissas(int timepoint) const {
		std::vector<double> abscissas;
		for (auto dim=0; dim<_domain.size(); dim++) {
			abscissas.push_back(_domain[dim]->get_ixn_param_traj()->get_ixn_param_at_timepoint(timepoint)->get_val());
		};
		return abscissas;
	};

	double DiffEqRHS::get_val_at_timepoint(int timepoint) const {
		return Grid::get_val(_form_abscissas(timepoint));
	};
	double DiffEqRHS::get_deriv_wrt_u_at_timepoint(int timepoint, q3c1::IdxSet global_vertex_idxs, std::vector<q3c1::DimType> dim_types) const {
		return Grid::get_deriv_wrt_coeff(_form_abscissas(timepoint), global_vertex_idxs, dim_types);
	};
	double DiffEqRHS::get_deriv_wrt_nu_at_timepoint(int timepoint, int deriv_dim) const {
		return Grid::get_deriv_wrt_abscissa(_form_abscissas(timepoint), deriv_dim);
	};

	/********************
	Update
	********************/

    void DiffEqRHS::update_calculate_and_store(int timepoint_start, int timepoint_end, double dt) {

		// Adjoint
		double adjoint_val;

		// Abscissas
		std::vector<double> abscissas;

		// Updates
		std::map<q3c1::Vertex*,std::vector<double>> updates;

		// Go through all times
		for (auto timepoint=timepoint_start; timepoint<timepoint_end; timepoint++) {

			// Adjoint
			adjoint_val = _parent_ixn_param_traj->get_adjoint()->get_val_at_timepoint(timepoint);

			// Form abscissas
			abscissas = _form_abscissas(timepoint);

			// Get updates for verts
			updates = Grid::get_deriv_wrt_coeffs_for_all_surrounding_verts(abscissas);

			// Append to any existing
			for (auto &pr: updates) {
				auto it = _updates.find(pr.first);
				if (it == _updates.end()) {
					// new
					_updates[pr.first] = pr.second;
					for (auto i=0; i<pr.second.size(); i++) {
						_updates[pr.first][i] *= dt * adjoint_val;
					};
				} else {
					// append
					for (auto i=0; i<pr.second.size(); i++) {
						_updates[pr.first][i] += dt * adjoint_val * pr.second[i];
					};
				};
			};
		};
	};

	// Committ the update
    void DiffEqRHS::update_committ_stored_sgd(double dopt) {
    };
    void DiffEqRHS::update_committ_stored_nesterov(double dopt, double nesterov_acc) {
    };
    void DiffEqRHS::update_committ_stored_adam(double dopt, int opt_step, double beta_1, double beta_2, double eps) {
        if (!_adam_m) {
            _adam_m = new std::map<q3c1::Vertex*,std::vector<double>>();
        };
        if (!_adam_v) {
            _adam_v = new std::map<q3c1::Vertex*,std::vector<double>>();
        };
        
        int opt_step_use = std::max(opt_step,1);
        
        double mhat, vhat;
        for (auto pr: _updates) {
            // _updates is never cleared
            // Therefore it will always contain the complete set of keys!
            // However, _adam_m and _adam_v may not contain the keys!
            
            // adam_m
            auto itm = _adam_m->find(pr.first);
            if (itm != _adam_m->end()) {
                // Exists, increment
                for (auto i=0; i<pr.second.size(); i++) {
                    (*_adam_m)[pr.first][i] = beta_1*itm->second[i] + (1.0 - beta_1)*pr.second[i];
                };
            } else {
                // Is zero
                (*_adam_m)[pr.first] = std::vector<double>(pr.second.size(),0.0);
                for (auto i=0; i<pr.second.size(); i++) {
                    (*_adam_m)[pr.first][i] = (1.0 - beta_1)*pr.second[i];
                };
            };
            
            // adam_v
            auto itv  = _adam_v->find(pr.first);
            if (itv != _adam_v->end()) {
                // Exists, increment
                for (auto i=0; i<pr.second.size(); i++) {
                    (*_adam_v)[pr.first][i] = beta_2*itm->second[i] + (1.0 - beta_2)*pow(pr.second[i],2);
                };
            } else {
                // Is zero
                (*_adam_v)[pr.first] = std::vector<double>(pr.second.size(),0.0);
                for (auto i=0; i<pr.second.size(); i++) {
                    (*_adam_v)[pr.first][i] = (1.0 - beta_2)*pow(pr.second[i],2);
                };
            };
            
            // Corrections and update
            for (auto i=0; i<pr.second.size(); i++) {
                // Corrections
                mhat = (*_adam_m)[pr.first][i] / (1.0 - pow(beta_1,opt_step_use));
                vhat = (*_adam_v)[pr.first][i] / (1.0 - pow(beta_2,opt_step_use));
                
                // Update val
                std::vector<q3c1::DimType> dims;
                if (_no_dims == 1) {
                    if (i==0) {
                        // Val
                        dims = {q3c1::DimType::VAL};
                    } else if (i==1) {
                        // Deriv
                        dims = {q3c1::DimType::DERIV};
                    };
                } else if (_no_dims == 2) {
                    if (i==0) {
                        // Val-val
                        dims = {q3c1::DimType::VAL,q3c1::DimType::VAL};
                    } else if (i==1) {
                        // Val-deriv
                        dims = {q3c1::DimType::VAL,q3c1::DimType::DERIV};
                    } else if (i==2) {
                        // Deriv-val
                        dims = {q3c1::DimType::DERIV,q3c1::DimType::VAL};
                    } else if (i==3) {
                        // Deriv-deriv
                        dims = {q3c1::DimType::DERIV,q3c1::DimType::DERIV};
                    };
                } else if (_no_dims == 3) {
                    if (i==0) {
                        // Val-val-val
                        dims = {q3c1::DimType::VAL,q3c1::DimType::VAL,q3c1::DimType::VAL};
                    } else if (i==1) {
                        // Val-val-deriv
                        dims = {q3c1::DimType::VAL,q3c1::DimType::VAL,q3c1::DimType::DERIV};
                    } else if (i==2) {
                        // Val-deriv-val
                        dims = {q3c1::DimType::VAL,q3c1::DimType::DERIV,q3c1::DimType::VAL};
                    } else if (i==3) {
                        // Deriv-val-val
                        dims = {q3c1::DimType::DERIV,q3c1::DimType::VAL,q3c1::DimType::VAL};
                    } else if (i==4) {
                        // Val-deriv-deriv
                        dims = {q3c1::DimType::VAL,q3c1::DimType::DERIV,q3c1::DimType::DERIV};
                    } else if (i==5) {
                        // Deriv-val-deriv
                        dims = {q3c1::DimType::DERIV,q3c1::DimType::VAL,q3c1::DimType::DERIV};
                    } else if (i==6) {
                        // Deriv-deriv-val
                        dims = {q3c1::DimType::DERIV,q3c1::DimType::DERIV,q3c1::DimType::VAL};
                    } else if (i==7) {
                        // Deriv-deriv-deriv
                        dims = {q3c1::DimType::DERIV,q3c1::DimType::DERIV,q3c1::DimType::DERIV};
                    };
                };
                
                // Update
                // _val -= dopt * mhat / (sqrt(vhat) + eps);
                pr.first->get_bf(dims)->increment_coeff(- dopt * mhat / (sqrt(vhat) + eps));
            };
            
            // Set the updates to zero
            // But do NOT clear the keys!
            _updates[pr.first] = std::vector<double>(pr.second.size(),0.0);
        };
    };
    
    // Verbose
    void DiffEqRHS::print_update_stored() const {
        for (auto const &pr: _updates) {
            std::cout << pr.first << ": " << std::flush;
            for (auto x: pr.second) {
                std::cout << x << " " << std::flush;
            };
            std::cout << std::endl;
        };
    };
};
