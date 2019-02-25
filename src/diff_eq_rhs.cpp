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
		_move(other);
	};
	Domain1D& Domain1D::operator=(const Domain1D& other) {
		if (this != &other)
		{
			_clean_up();
            Dimension1D::operator=(other);
            _copy(other);
		};
		return *this;
	};
	Domain1D& Domain1D::operator=(Domain1D&& other) {
		if (this != &other)
		{
			_clean_up();
            Dimension1D::operator=(std::move(other));
            _move(other);
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
	void Domain1D::_move(Domain1D& other)
	{
        _ixn_param_traj = other._ixn_param_traj;
        other._ixn_param_traj = nullptr;
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

	DiffEqRHS::DiffEqRHS(std::string name, ITptr parent_ixn_param_traj, std::vector<Domain1D*> domain, double lr) : q3c1::Grid(std::vector<q3c1::Dimension1D*>(domain.begin(),domain.end())) {
		_name = name;
		_domain = domain;
		_no_dims = _domain.size();
		_parent_ixn_param_traj = parent_ixn_param_traj;
        _abscissas = std::vector<double>(_no_dims,0.0);
        _lr = lr;
        
		if (_no_dims == 0 || _no_dims > 3) {
			std::cerr << ">>> Error: DiffEqRHS::DiffEqRHS <<< Only dims 1,2,3 are supported" << std::endl;
			exit(EXIT_FAILURE);
		};
        
        // Number coeffs
        _no_coeffs = pow(2,_no_dims);
        
        // Coefficient order
        if (_no_dims == 1) {
            // Val
            _coeff_order.push_back({q3c1::DimType::VAL});
            // Deriv
            _coeff_order.push_back({q3c1::DimType::DERIV});
        } else if (_no_dims == 2) {
            // Val-val
            _coeff_order.push_back({q3c1::DimType::VAL,q3c1::DimType::VAL});
            // Val-Deriv
            _coeff_order.push_back({q3c1::DimType::VAL,q3c1::DimType::DERIV});
            // Deriv-Val
            _coeff_order.push_back({q3c1::DimType::DERIV,q3c1::DimType::VAL});
            // Deriv-Deriv
            _coeff_order.push_back({q3c1::DimType::DERIV,q3c1::DimType::DERIV});
        } else if (_no_dims == 3) {
            // Val-val-val
            _coeff_order.push_back({q3c1::DimType::VAL,q3c1::DimType::VAL,q3c1::DimType::VAL});
            // Val-Val-Deriv
            _coeff_order.push_back({q3c1::DimType::VAL,q3c1::DimType::VAL,q3c1::DimType::DERIV});
            // Val-Deriv-Val
            _coeff_order.push_back({q3c1::DimType::VAL,q3c1::DimType::DERIV,q3c1::DimType::VAL});
            // Deriv-Val-Val
            _coeff_order.push_back({q3c1::DimType::DERIV,q3c1::DimType::VAL,q3c1::DimType::VAL});
            // Val-deriv-deriv
            _coeff_order.push_back({q3c1::DimType::VAL,q3c1::DimType::DERIV,q3c1::DimType::DERIV});
            // Deriv-Val-Deriv
            _coeff_order.push_back({q3c1::DimType::DERIV,q3c1::DimType::VAL,q3c1::DimType::DERIV});
            // Deriv-Deriv-Val
            _coeff_order.push_back({q3c1::DimType::DERIV,q3c1::DimType::DERIV,q3c1::DimType::VAL});
            // Deriv-Deriv-Deriv
            _coeff_order.push_back({q3c1::DimType::DERIV,q3c1::DimType::DERIV,q3c1::DimType::DERIV});
        };
        
        // Nullptr for adams/nesterovs
        _nesterov_y_s = nullptr;
        _nesterov_y_sp1 = nullptr;
        _adam_m = nullptr;
        _adam_v = nullptr;
        
        // Max update
        _mag_max_update = nullptr;
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
        _no_coeffs = other._no_coeffs;
        _coeff_order = other._coeff_order;
		_domain = other._domain;
		_parent_ixn_param_traj = other._parent_ixn_param_traj;
		_updates = other._updates;
        _abscissas = other._abscissas;
        _lr = other._lr;
        
        if (other._nesterov_y_s) {
			_nesterov_y_s = new std::map<q3c1::Vertex*,std::vector<double>>(*other._nesterov_y_s);
        } else {
            _nesterov_y_s = nullptr;
        };
		if (other._nesterov_y_sp1) {
			_nesterov_y_sp1 = new std::map<q3c1::Vertex*,std::vector<double>>(*other._nesterov_y_sp1);
        } else {
            _nesterov_y_sp1 = nullptr;
        };
        if (other._adam_m) {
            _adam_m = new std::map<q3c1::Vertex*,std::vector<double>>(*other._adam_m);
        } else {
            _adam_m = nullptr;
        };
        if (other._adam_v) {
            _adam_v = new std::map<q3c1::Vertex*,std::vector<double>>(*other._adam_v);
        } else {
            _adam_v = nullptr;
        };
        
        if (other._mag_max_update) {
            _mag_max_update = new double(*other._mag_max_update);
        } else {
            _mag_max_update = nullptr;
        };
	};
	void DiffEqRHS::_move(DiffEqRHS& other) {
		_name = other._name;
		_no_dims = other._no_dims;
        _no_coeffs = other._no_coeffs;
        _coeff_order = other._coeff_order;
		_domain = other._domain;
		_parent_ixn_param_traj = other._parent_ixn_param_traj;
		_updates = other._updates;
        _abscissas = other._abscissas;
        _lr = other._lr;

		_nesterov_y_s = other._nesterov_y_s;
		_nesterov_y_sp1 = other._nesterov_y_sp1;

        _adam_m = other._adam_m;
        _adam_v = other._adam_v;
        
        _mag_max_update = other._mag_max_update;
        
		// Reset other
		other._name = "";
		other._no_dims = 0;
        other._no_coeffs = 0;
        other._coeff_order.clear();
		other._domain.clear();
		other._parent_ixn_param_traj = nullptr;
		other._updates.clear();
        other._abscissas.clear();
        other._lr = 0.0;
        
        other._nesterov_y_s = nullptr;
		other._nesterov_y_sp1 = nullptr;
        
        other._adam_m = nullptr;
        other._adam_v = nullptr;
        
        other._mag_max_update = nullptr;
	};

	void DiffEqRHS::_clean_up() {
		// Note: Domain is not owned => do not delete
		if (_nesterov_y_s) {
			delete _nesterov_y_s;
		};
        _nesterov_y_s = nullptr;
        
        if (_nesterov_y_sp1) {
			delete _nesterov_y_sp1;
		};
        _nesterov_y_sp1 = nullptr;
        
        if (_adam_m) {
            delete _adam_m;
        };
        _adam_m = nullptr;
        
        if (_adam_v) {
            delete _adam_v;
        };
        _adam_v = nullptr;
        
        if (_mag_max_update) {
            delete _mag_max_update;
        };
        _mag_max_update = nullptr;
	};

    // ***************
    // MARK: - Learning rate
    // ***************
    
    void DiffEqRHS::set_lr(double lr) {
        _lr = lr;
    };
    double DiffEqRHS::get_lr() const {
        return _lr;
    };

	/********************
	Validate
	********************/

	// Validate setup
	void DiffEqRHS::validate_setup() const {
		std::cout << "--- Validate Basis func " << _name << " ---" << std::endl;
	};

    // ***************
    // MARK: - Set magnitude of max update
    // ***************
    
    void DiffEqRHS::set_mag_max_update(double mag) {
        if (_mag_max_update) {
            delete _mag_max_update;
        };
        
        if (mag < 0.0) {
            std::cerr << ">>> DiffEqRHS::set_mag_max_update <<< Error: mag of max update must be positive" << std::endl;
            exit(EXIT_FAILURE);
        };
        
        _mag_max_update = new double(mag);
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

    void DiffEqRHS::_form_abscissas(int timepoint) const {
		for (auto dim=0; dim<_no_dims; dim++) {
            _abscissas[dim] = _domain[dim]->get_ixn_param_traj()->get_ixn_param_at_timepoint(timepoint)->get_val();
		};
	};
    void DiffEqRHS::_form_abscissas(const std::map<ITptr,std::vector<double>>& vals, int idx) const {
        for (auto dim=0; dim<_no_dims; dim++) {
            _abscissas[dim] = vals.at(_domain[dim]->get_ixn_param_traj()).at(idx);
        };
    };

	double DiffEqRHS::get_val_at_timepoint(int timepoint) const {
        _form_abscissas(timepoint);
        return Grid::get_val(_abscissas);
	};
    double DiffEqRHS::get_val_from_map(const std::map<ITptr,std::vector<double>>& vals, int idx) const {
        _form_abscissas(vals, idx);
        return Grid::get_val(_abscissas);
    };
    
	double DiffEqRHS::get_deriv_wrt_u_at_timepoint(int timepoint, q3c1::IdxSet global_vertex_idxs, std::vector<q3c1::DimType> dim_types) const {
        _form_abscissas(timepoint);
		return Grid::get_deriv_wrt_coeff(_abscissas, global_vertex_idxs, dim_types);
	};
	double DiffEqRHS::get_deriv_wrt_nu_at_timepoint(int timepoint, int deriv_dim) const {
        _form_abscissas(timepoint);
		return Grid::get_deriv_wrt_abscissa(_abscissas, deriv_dim);
	};

    // ***************
    // MARK: - Fix vertices at some timepoint
    // ***************
    
    void DiffEqRHS::fix_all_verts_around_at_timepoint(int timepoint, bool fixed) const {
        _form_abscissas(timepoint);
        std::pair<q3c1::Cell*,std::vector<double>> cell_pr = get_cell(_abscissas);
        for (auto vert_pr: cell_pr.first->get_all_vertices()) {
            for (auto bf: vert_pr.second->get_bfs()) {
                bf->set_is_val_fixed(fixed);
            };
        };
    };
    
	/********************
	Update
	********************/

    void DiffEqRHS::update_calculate_and_store(int timepoint_start, int timepoint_end, double dt) {

        // Clear before starting
        _updates.clear();
        
		// Adjoint
		double adjoint_val;

		// Updates
		std::map<q3c1::Vertex*,std::vector<double>> updates;
        
		// Go through all times
		for (auto timepoint=timepoint_start; timepoint<=timepoint_end; timepoint++) {

			// Adjoint
			adjoint_val = _parent_ixn_param_traj->get_adjoint()->get_val_at_timepoint(timepoint);

			// Form abscissas
            _form_abscissas(timepoint);

			// Get updates for verts
			updates = Grid::get_deriv_wrt_coeffs_for_all_surrounding_verts(_abscissas);

			// Append to any existing
			for (auto &pr: updates) {
				auto it = _updates.find(pr.first);
				if (it == _updates.end()) {
					// new
                    _updates[pr.first] = std::vector<double>(_no_coeffs,0.0);
                };
                // append
                for (auto i=0; i<_no_coeffs; i++) {
                    _updates[pr.first][i] -= dt * adjoint_val * pr.second[i];
                    
                    // std::cout << "update_calculate_and_store: timepoint = " << timepoint << " key ptr: " << pr.first << " idx: " << i << " val: " << dt * adjoint_val * pr.second[i] << std::endl;
                    
                };
			};
		};
	};

	// Committ the update
    void DiffEqRHS::update_committ_stored_sgd() {
        if (_mag_max_update) {
            
            for (auto pr: _updates) {
                for (auto i=0; i<_no_coeffs; i++) {
                    if (-_lr * pr.second.at(i) > (*_mag_max_update)) {
                        pr.first->get_bf(_coeff_order.at(i))->increment_coeff(*_mag_max_update);
                    } else if (-_lr * pr.second.at(i) < -1.0 * (*_mag_max_update)) {
                        pr.first->get_bf(_coeff_order.at(i))->increment_coeff(-1.0 * *_mag_max_update);
                    } else {
                        pr.first->get_bf(_coeff_order.at(i))->increment_coeff(- _lr * pr.second.at(i));
                    };
                };
            };
            
        } else {
            
            // double ave_mag=0.0,max_mag=0.0;
            // int counts = 0;
            for (auto pr: _updates) {
                for (auto i=0; i<_no_coeffs; i++) {
                    pr.first->get_bf(_coeff_order.at(i))->increment_coeff(- _lr * pr.second.at(i));
                    /*
                    if (i==0) {
                        ave_mag += abs( _lr * pr.second.at(i) );
                        max_mag = std::max(max_mag,abs( _lr * pr.second.at(i) ));
                        counts++;
                    };
                     */
                };
            };
            // ave_mag /= counts;
            
            // std::cout << _name << ": average magnitude of updates: " << ave_mag << std::endl;
            // std::cout << _name << ": max magnitude of updates: " << max_mag << std::endl;
        };
    };
    void DiffEqRHS::update_committ_stored_nesterov(double nesterov_acc) {
        std::cerr << ">>> DiffEqRHS::update_committ_stored_nesterov <<< Unsupported currently" << std::endl;
        exit(EXIT_FAILURE);
    };
    void DiffEqRHS::update_committ_stored_adam(int opt_step, double beta_1, double beta_2, double eps) {

        // Create adam variables if necessary
        if (!_adam_m) {
            _adam_m = new std::map<q3c1::Vertex*,std::vector<double>>();
            // std::cout << "created _adam_m = " << _adam_m << std::endl;
        };
        if (!_adam_v) {
            _adam_v = new std::map<q3c1::Vertex*,std::vector<double>>();
        };
        
        int opt_step_use = std::max(opt_step,1);
        
        double mhat, vhat, update_val;
        std::map<q3c1::Vertex*,std::vector<double>>::iterator itm, itv;
        for (auto pr: _updates) {
            
            // adam_m
            itm = _adam_m->find(pr.first);
            if (itm != _adam_m->end()) {
                // Exists, increment
                for (auto i=0; i<_no_coeffs; i++) {
                    (*_adam_m)[pr.first][i] = beta_1*itm->second[i] + (1.0 - beta_1)*pr.second[i];
                };
            } else {
                // Is zero
                (*_adam_m)[pr.first] = std::vector<double>(_no_coeffs,0.0);
                for (auto i=0; i<_no_coeffs; i++) {
                    (*_adam_m)[pr.first][i] = (1.0 - beta_1)*pr.second[i];
                };
            };
            
            // adam_v
            itv  = _adam_v->find(pr.first);
            if (itv != _adam_v->end()) {
                // Exists, increment
                for (auto i=0; i<_no_coeffs; i++) {
                    (*_adam_v)[pr.first][i] = beta_2*itv->second[i] + (1.0 - beta_2)*pow(pr.second[i],2);
                };
            } else {
                // Is zero
                (*_adam_v)[pr.first] = std::vector<double>(_no_coeffs,0.0);
                for (auto i=0; i<_no_coeffs; i++) {
                    (*_adam_v)[pr.first][i] = (1.0 - beta_2)*pow(pr.second[i],2);
                };
            };
    
            /*
            std::cout << "adams are:" << std::endl;
            for (auto i=0; i<_no_coeffs; i++) {
                std::cout << (*_adam_m)[pr.first][i] << " ";
            };
            std::cout << std::endl;
             */
            
            // Corrections and update
            for (auto i=0; i<_no_coeffs; i++) {
                
                // Corrections
                mhat = (*_adam_m)[pr.first][i] / (1.0 - pow(beta_1,opt_step_use));
                vhat = (*_adam_v)[pr.first][i] / (1.0 - pow(beta_2,opt_step_use));
                
                // Update
                // _val -= dopt * mhat / (sqrt(vhat) + eps);
                update_val = - _lr * mhat / (sqrt(vhat) + eps);
                if (_mag_max_update) {
                    if (update_val > *_mag_max_update) {
                        std::cout << "Warning: update: " << update_val << " is beyond the max magnitude; clipping to: " << *_mag_max_update << std::endl;
                        update_val = *_mag_max_update;
                    } else if (update_val < - *_mag_max_update) {
                        std::cout << "Warning: update: " << update_val << " is beyond the max magnitude; clipping to: " << - *_mag_max_update << std::endl;
                        update_val = - *_mag_max_update;
                    };
                };
                
                pr.first->get_bf(_coeff_order.at(i))->increment_coeff(update_val);
                // std::cout << "update_committ_stored_adam: ptr: " << pr.first << " idx: " << i << " update: " << - dopt * mhat / (sqrt(vhat) + eps) << std::endl;
            };
            
            // Set the updates to zero
            // But do NOT clear the keys!
            // _updates[pr.first] = std::vector<double>(_no_coeffs,0.0);
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
