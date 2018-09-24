#include <dynamicboltz>

#include <iostream>

using namespace dblz;
using namespace std;
using namespace dcu;

int main() {

	/****************************************
	Make some species
	****************************************/

	cout << "--- Making species ---" << endl;

	auto species_A = make_shared<Species>("A");

	cout << "--- [Finished] Making species ---" << endl;
	cout << endl;

	/****************************************
	Make ixn
	****************************************/

	cout << "--- Making ixn func ---" << endl;

	// J init = 0.1
	// h init = 0.5
	double init_val = 0.5;
	auto ixn = make_shared<IxnParam>("bias for A",IxnParamType::H, init_val);

	cout << "--- [Finished] Making ixn func ---" << endl;
	cout << endl;

	/****************************************
	Make RHS diff eq for ixn
	****************************************/

	cout << "--- Making RHS diff eq ---" << endl;

	// Domain
	double min_val=-1.0, max_val=0.6;
	int no_pts=21;
	Domain1D domain_1d = Domain1D(ixn,min_val,max_val,no_pts);
	Domain domain({domain_1d});

	// RHS
	auto rhs = make_shared<DiffEqRHS>("bias DE RHS", ixn, domain);

	// Add to ixn
	ixn->set_diff_eq_rhs(rhs);

	cout << "--- [Finished] Making RHS diff eq ---" << endl;
	cout << endl;

	/****************************************
	Adjoint
	****************************************/

	cout << "--- Making adjoint ---" << endl;

	auto adjoint = make_shared<Adjoint>("adjoint",ixn);
	ixn->set_adjoint(adjoint);

	cout << "--- [Finished] Making adjoint ---" << endl;
	cout << endl;

	/****************************************
	Setup lattice
	****************************************/

	cout << "--- Making lattice ---" << endl;

	// 10x10x10 cube
	auto latt = make_shared<Lattice>(3,10);

	// Set possible species
	latt->all_units_v_add_possible_species(species_A);

	// Make ixn dict
	auto bias_dict = make_shared<BiasDict>();
	bias_dict->add_ixn(species_A,ixn);

	// Add ixn to lattice
	latt->all_units_v_set_bias_dict(bias_dict);

	// Add all units to moment of ixn func
	for (auto i1=1; i1<= 10; i1++) {
		for (auto i2=1; i2<= 10; i2++) {
			for (auto i3=1; i3<= 10; i3++) {
				ixn->get_moment()->add_unit_to_monitor_h(&latt->get_unit_v(i1,i2,i3));
			};
		};
	};

	cout << "--- [Finished] Making lattice ---" << endl;
	cout << endl;

	/****************************************
	Do an iteration of the learning problem
	****************************************/

	// Opt solver class
	OptProblem opt(latt);

	// Add ixn param
	opt.add_ixn_param(ixn);

	// Params
	int no_opt_steps = 50;
	int no_timesteps = 20;
	double dt=0.1;
	int no_latt_sampling_steps = 10;
	int batch_size = 5;
	double dopt=0.002;

	// Filenames
	FNameSeriesColl fnames;
	for (auto i_batch=0; i_batch<batch_size; i_batch++) {

		// Series of filenames
		FNameSeries fname_series;
		for (auto timepoint=0; timepoint<=no_timesteps; timepoint++) {
			fname_series.fnames.push_back("stoch_sim/lattice_v"+pad_str(i_batch,2)+"/lattice/"+pad_str(timepoint,4)+".txt");
		};

		// Add
		fnames.add_fname_series(fname_series);
	};

	// Options
	OptionsSolve options;
	options.VERBOSE_NU = false;
	options.VERBOSE_WAKE_ASLEEP = false;


	// Solve
	opt.solve(no_opt_steps,no_timesteps,batch_size,dt,dopt,no_latt_sampling_steps,fnames,options);

	return 0;
};