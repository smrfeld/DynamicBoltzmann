#include <dynamicboltz>

#include <iostream>

using namespace dblz;
using namespace std;
using namespace dcu;

int main() {

	/****************************************
	Make some species
	****************************************/

	auto species_A = make_shared<Species>("A");

	cout << "Made species: " << species_A->get_name() << endl;
	cout << endl;

	/****************************************
	Make ixn
	****************************************/

	double init_val = 0.83;
	auto ixn = make_shared<IxnParam>("bias for A",IxnParamType::H, init_val);

	cout << "Made ixn bias for A: " << ixn->get_name() << endl;
	cout << endl;

	/****************************************
	Make RHS diff eq for ixn
	****************************************/

	// Domain
	double min_val=-1.0, max_val=1.0;
	int no_pts=21;
	Domain1D domain_1d = Domain1D(ixn,min_val,max_val,no_pts);
	Domain domain({domain_1d});

	// RHS
	auto rhs = make_shared<DiffEqRHS>("bias DE RHS", ixn, domain);

	// Add to ixn
	ixn->set_diff_eq_rhs(rhs);

	cout << "Made RHS of diff eq: " << rhs->get_name() << endl;
	cout << endl;

	/****************************************
	Set some random grid points
	****************************************/

	std::vector<int> idx_set({0});
	idx_set[0] = 15;
	rhs->get_grid_point_ref(idx_set).set_ordinate(-0.1);
	idx_set[0] = 16;
	rhs->get_grid_point_ref(idx_set).set_ordinate(-0.2);
	idx_set[0] = 17;
	rhs->get_grid_point_ref(idx_set).set_ordinate(-0.3);
	idx_set[0] = 18;
	rhs->get_grid_point_ref(idx_set).set_ordinate(-0.4);
	idx_set[0] = 19;
	rhs->get_grid_point_ref(idx_set).set_ordinate(-0.5);
	idx_set[0] = 20;
	rhs->get_grid_point_ref(idx_set).set_ordinate(-0.6);

	// rhs->write_vals("test_rhs_diff_eq.txt");

	cout << "Set & wrote RHS" << endl;
	cout << endl;

	/****************************************
	Get value of diff eq
	****************************************/

	std::cout << "Val at time 0 = " << rhs->get_val_at_timepoint(0) << std::endl;
	std::cout << "Deriv wrt nu at time 0 = " << rhs->get_deriv_wrt_nu_at_timepoint(0,0) << std::endl;

	/****************************************
	Solve
	****************************************/

	ixn->set_no_timesteps(20);
	double dt=0.1;

	std::cout << "t = " << 0 << " val = " << ixn->get_val_at_timepoint(0) << std::endl;
	for (auto t=0; t<20; t++) {
		ixn->solve_diff_eq_at_timepoint_to_plus_one(t,dt);
		std::cout << "t = " << t+1 << " val = " << ixn->get_val_at_timepoint(t+1) << std::endl;
	};

	cout << "Solved diff eq for state var" << endl;
	cout << endl;

	/****************************************
	Adjoint
	****************************************/

	auto adjoint = make_shared<Adjoint>("adjoint",ixn);
	ixn->set_adjoint(adjoint);

	std::cout << "Made adjoint: " << adjoint->get_name() << std::endl;
	std::cout << std::endl;

	/****************************************
	Solve
	****************************************/

	adjoint->set_no_timesteps(20);
	adjoint->set_zero_end_cond_timepoint(20);
	std::cout << "t = " << 20 << " val = " << adjoint->get_val_at_timepoint(20) << std::endl;
	for (auto t=20; t>=1; t--) {
		adjoint->solve_diff_eq_at_timepoint_to_minus_one(t,dt);
		std::cout << "t = " << t-1 << " val = " << adjoint->get_val_at_timepoint(t-1) << std::endl;
	};

	cout << "Solved diff eq for adjoint" << endl;
	cout << endl;

	/****************************************
	Form the update
	****************************************/

	double dopt=0.1;
	rhs->calculate_update(0,20,dt,dopt);

	cout << "Formed update" << endl;
	rhs->print_update_stored();
	cout << endl;

	return 0;
};




