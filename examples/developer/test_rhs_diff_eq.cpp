#include <dynamicboltz>

#include <iostream>

using namespace dblz;
using namespace std;

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
	vector<Domain1D> domain({Domain1D(ixn,min_val,max_val,no_pts)});

	// RHS
	auto rhs = make_shared<DiffEqRHS>("bias DE RHS",domain);

	// Add to ixn
	ixn->set_diff_eq(rhs);

	cout << "Made RHS of diff eq: " << rhs->get_name() << endl;
	cout << endl;

	/****************************************
	Set some random grid points
	****************************************/

	rhs->set_by_idx(15,-0.1);
	rhs->set_by_idx(16,-0.2);
	rhs->set_by_idx(17,-0.3);
	rhs->set_by_idx(18,-0.4);
	rhs->set_by_idx(19,-0.5);
	rhs->set_by_idx(20,-0.6);

	rhs->write_vals("test_rhs_diff_eq.txt");

	cout << "Set & wrote RHS" << endl;
	cout << endl;

	/****************************************
	Get value of diff eq
	****************************************/

	std::cout << "Val at time 0 = " << rhs->get_val_at_timepoint(0) << std::endl;
	std::cout << "Deriv at time 0 = " << rhs->get_deriv_at_timepoint(0,0) << std::endl;

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

	cout << "Solved diff eq" << endl;
	cout << endl;

	return 0;
};




