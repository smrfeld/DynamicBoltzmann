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
	auto species_B = make_shared<Species>("B");

	cout << "--- [Finished] Making species ---" << endl;
	cout << endl;

	/****************************************
	Make ixn
	****************************************/

	cout << "--- Making ixn func ---" << endl;

	auto ixn_hA = make_shared<IxnParam>("hA",IxnParamType::H, 1.0);
	auto ixn_hB = make_shared<IxnParam>("hB",IxnParamType::H, 1.0);
	auto ixn_jAB = make_shared<IxnParam>("jAB",IxnParamType::J, 0.0);

	cout << "--- [Finished] Making ixn func ---" << endl;
	cout << endl;

	/****************************************
	Make RHS diff eq for ixn
	****************************************/

	cout << "--- Making RHS diff eq ---" << endl;

	// Domain
	Domain1D domain_1d_hA = Domain1D(ixn_hA,-3.5,1.2,16);
	Domain1D domain_1d_hB = Domain1D(ixn_hB,-3.5,1.2,16);
	Domain domain({domain_1d_hA,domain_1d_hB});

	// RHS
	auto rhs_hA = make_shared<DiffEqRHS>("DE hA", ixn_hA, domain);
	ixn_hA->set_diff_eq_rhs(rhs_hA);
	auto rhs_hB = make_shared<DiffEqRHS>("DE hB", ixn_hB, domain);
	ixn_hB->set_diff_eq_rhs(rhs_hB);
	auto rhs_jAB = make_shared<DiffEqRHS>("DE jAB", ixn_jAB, domain);
	ixn_jAB->set_diff_eq_rhs(rhs_jAB);

	cout << "--- [Finished] Making RHS diff eq ---" << endl;
	cout << endl;

	/****************************************
	Adjoint
	****************************************/

	cout << "--- Making adjoint ---" << endl;

	auto adjoint_hA = make_shared<Adjoint>("adjoint hA",ixn_hA);
	ixn_hA->set_adjoint(adjoint_hA);

	auto adjoint_hB = make_shared<Adjoint>("adjoint hB",ixn_hB);
	ixn_hB->set_adjoint(adjoint_hB);

	auto adjoint_jAB = make_shared<Adjoint>("adjoint jAB",ixn_jAB);
	ixn_jAB->set_adjoint(adjoint_jAB);

	cout << "--- [Finished] Making adjoint ---" << endl;
	cout << endl;

	/****************************************
	Setup lattice
	****************************************/

	cout << "--- Making lattice ---" << endl;

	// 10x10x10 cube
	auto latt = make_shared<Lattice>(3,10);

	// NN connectivity
	latt->all_conns_vv_init();

	// Set possible species
	latt->all_units_v_add_possible_species(species_A);
	latt->all_units_v_add_possible_species(species_B);

	// Make ixn dicts
	auto bias_dict = make_shared<BiasDict>();
	bias_dict->add_ixn(species_A,ixn_hA);
	bias_dict->add_ixn(species_B,ixn_hB);

	auto ixn_dict = make_shared<O2IxnDict>();
	ixn_dict->add_ixn(species_A,species_B,ixn_jAB);

	// Add ixn to lattice
	latt->all_units_v_set_bias_dict(bias_dict);
	latt->all_conns_vv_set_ixn_dict(ixn_dict);

	// Add all units to moment of ixn func
	latt->all_units_v_add_to_moment_h(ixn_hA->get_moment());
	latt->all_units_v_add_to_moment_h(ixn_hB->get_moment());
	latt->all_conns_vv_add_to_moment_j(ixn_jAB->get_moment());

	cout << "--- [Finished] Making lattice ---" << endl;
	cout << endl;

	/****************************************
	Do an iteration of the learning problem
	****************************************/

	// Opt solver class
	OptProblem opt(latt);

	// Add ixn params
	opt.add_ixn_param(ixn_hA);
	opt.add_ixn_param(ixn_hB);
	opt.add_ixn_param(ixn_jAB);

	// Params
	int no_opt_steps = 100;
	int no_timesteps = 100;
	double dt=0.1;
	int no_latt_sampling_steps = 10;
	int batch_size = 5;
	double dopt=0.002;

	// Filenames
	FNameSeriesColl fnames;
	for (auto i_batch=0; i_batch<10; i_batch++) {

		// Series of filenames
		FNameSeries fname_series;
		for (auto timepoint=0; timepoint<=no_timesteps; timepoint++) {
			fname_series.fnames.push_back("stoch_sim/lattice_v"+pad_str(i_batch,3)+"/lattice/"+pad_str(timepoint,4)+".txt");
		};

		// Add
		fnames.add_fname_series(fname_series);
	};

	// Options
	OptionsSolve options;
	options.VERBOSE_NU = false;
	options.VERBOSE_WAKE_ASLEEP = false;
	options.VERBOSE_MOMENT = false;
	options.MODE_random_integral_range = true;
	options.VAL_random_integral_range_size = 10;

	// Solve
	// opt.solve(no_opt_steps,no_timesteps,batch_size,dt,dopt,no_latt_sampling_steps,fnames,options);

	// Manually
	opt.init_structures(no_timesteps,batch_size);
	for (auto opt_step=1; opt_step<=no_opt_steps; opt_step++) {

		std::cout << "------------------" << std::endl;
		std::cout << "Opt step: " << opt_step << " / " << no_opt_steps << std::endl;
		std::cout << "------------------" << std::endl;

		// Take a step
		opt.solve_one_step(no_timesteps,batch_size,dt,dopt,no_latt_sampling_steps,fnames,options);

		if (opt_step <= 50) {
			options.VAL_random_integral_range_size = 10;
		} else if (opt_step <= 75) {
			options.VAL_random_integral_range_size = 25;
		} else {
			options.VAL_random_integral_range_size = 50;
		};

		// Sample moments (whole trajectory)
		if (opt_step==no_opt_steps) {
			opt.wake_sleep_loop(0,no_timesteps,batch_size,no_latt_sampling_steps,fnames,false);
		};

		// Print moments
		std::cout << ixn_hA->get_name() << " " << std::flush;
		ixn_hA->get_moment()->print_moment_comparison();

		std::cout << ixn_hB->get_name() << " " << std::flush;
		ixn_hB->get_moment()->print_moment_comparison();

		std::cout << ixn_jAB->get_name() << " " << std::flush;
		ixn_jAB->get_moment()->print_moment_comparison();

		// Write moments
		ixn_hA->get_moment()->write_to_file("data_learned/moments/"+ixn_hA->get_name()+"_"+pad_str(opt_step,3)+".txt");
		ixn_hB->get_moment()->write_to_file("data_learned/moments/"+ixn_hB->get_name()+"_"+pad_str(opt_step,3)+".txt");
		ixn_jAB->get_moment()->write_to_file("data_learned/moments/"+ixn_jAB->get_name()+"_"+pad_str(opt_step,3)+".txt");

		// Write ixn params
		ixn_hA->write_to_file("data_learned/ixn_params/"+ixn_hA->get_name()+"_"+pad_str(opt_step,3)+".txt");
		ixn_hB->write_to_file("data_learned/ixn_params/"+ixn_hB->get_name()+"_"+pad_str(opt_step,3)+".txt");
		ixn_jAB->write_to_file("data_learned/ixn_params/"+ixn_jAB->get_name()+"_"+pad_str(opt_step,3)+".txt");

		// Write adjoint
		ixn_hA->get_adjoint()->write_to_file("data_learned/adjoint/"+ixn_hA->get_name()+"_"+pad_str(opt_step,3)+".txt");
		ixn_hB->get_adjoint()->write_to_file("data_learned/adjoint/"+ixn_hB->get_name()+"_"+pad_str(opt_step,3)+".txt");
		ixn_jAB->get_adjoint()->write_to_file("data_learned/adjoint/"+ixn_jAB->get_name()+"_"+pad_str(opt_step,3)+".txt");

		// Write diff eq rhs
		if (opt_step == no_opt_steps) {
			rhs_hA->write_to_file("data_learned/diff_eq_rhs/"+ixn_hA->get_name()+"_"+pad_str(opt_step,3)+".txt");
			rhs_hB->write_to_file("data_learned/diff_eq_rhs/"+ixn_hB->get_name()+"_"+pad_str(opt_step,3)+".txt");
			rhs_jAB->write_to_file("data_learned/diff_eq_rhs/"+ixn_jAB->get_name()+"_"+pad_str(opt_step,3)+".txt");
		};
	};

	return 0;
};