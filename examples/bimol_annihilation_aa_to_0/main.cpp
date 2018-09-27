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
	auto species_X = make_shared<Species>("X");

	cout << "--- [Finished] Making species ---" << endl;
	cout << endl;

	/****************************************
	Make ixn
	****************************************/

	cout << "--- Making ixn func ---" << endl;

	auto ixn_hA = make_shared<IxnParam>("hA",IxnParamType::H, 1.0);
	auto ixn_wAX = make_shared<IxnParam>("wAX",IxnParamType::W, 0.0);
	auto ixn_bX = make_shared<IxnParam>("bX",IxnParamType::B, 0.0);

	cout << "--- [Finished] Making ixn func ---" << endl;
	cout << endl;

	/****************************************
	Make RHS diff eq for ixn
	****************************************/

	cout << "--- Making RHS diff eq ---" << endl;

	// Domain
	Domain1D domain_1d_hA = Domain1D(ixn_hA,-2.0,2.0,11);
	Domain1D domain_1d_wAX = Domain1D(ixn_wAX,-5.0,0.5,11);
	Domain1D domain_1d_bX = Domain1D(ixn_bX,-0.5,2.0,11);
	Domain domain({domain_1d_hA,domain_1d_wAX,domain_1d_bX});

	// RHS
	auto rhs_hA = make_shared<DiffEqRHS>("DE hA", ixn_hA, domain);
	ixn_hA->set_diff_eq_rhs(rhs_hA);
	auto rhs_wAX = make_shared<DiffEqRHS>("DE wAX", ixn_wAX, domain);
	ixn_wAX->set_diff_eq_rhs(rhs_wAX);
	auto rhs_bX = make_shared<DiffEqRHS>("DE bX", ixn_bX, domain);
	ixn_bX->set_diff_eq_rhs(rhs_bX);

	cout << "--- [Finished] Making RHS diff eq ---" << endl;
	cout << endl;

	/****************************************
	Adjoint
	****************************************/

	cout << "--- Making adjoint ---" << endl;

	auto adjoint_hA = make_shared<Adjoint>("adjoint hA",ixn_hA);
	ixn_hA->set_adjoint(adjoint_hA);

	auto adjoint_wAX = make_shared<Adjoint>("adjoint wAX",ixn_wAX);
	ixn_wAX->set_adjoint(adjoint_wAX);

	auto adjoint_bX = make_shared<Adjoint>("adjoint bX",ixn_bX);
	ixn_bX->set_adjoint(adjoint_bX);

	cout << "--- [Finished] Making adjoint ---" << endl;
	cout << endl;

	/****************************************
	Setup lattice
	****************************************/

	cout << "--- Making lattice ---" << endl;

	// 1000 cube
	auto latt = make_shared<Lattice>(1,1000);

	/*****
	Visible
	*****/

	cout << " > begin visible" << endl;

	// Set possible species
	latt->all_units_v_add_possible_species(species_A);

	// Make bias dict
	auto bias_dict_visible = make_shared<BiasDict>();
	bias_dict_visible->add_ixn(species_A,ixn_hA);

	// Add dict to latt
	latt->all_units_v_set_bias_dict(bias_dict_visible);

	// Add all units to moment of ixn func
	latt->all_units_v_add_to_moment_h(ixn_hA->get_moment());

	cout << " > end visible" << endl;

	/*****
	Hidden
	*****/

	cout << " > begin hidden" << endl;

	// Add hidden units
	int layer = 1;
	for (auto i=1; i<1000; i++) { // 999 total
		latt->add_hidden_unit(layer,i);
	};

	// Add species
	latt->all_units_h_add_possible_species(species_X);

	// Connect hidden units
	for (auto i=1; i<1000; i++) {
		// Connects i with i+1
		latt->add_conn_vh(&latt->get_unit_v(i),&latt->get_unit_h(layer,i));
		latt->add_conn_vh(&latt->get_unit_v(i+1),&latt->get_unit_h(layer,i));
	};

	// Make bias dict
	auto bias_dict_hidden = make_shared<BiasDict>();
	bias_dict_hidden->add_ixn(species_X,ixn_bX);

	// Add dict to latt
	latt->all_units_h_set_bias_dict(bias_dict_hidden);

	// Make ixn dict for W
	auto ixn_dict = make_shared<O2IxnDict>();
	ixn_dict->add_ixn(species_A,species_X,ixn_wAX);

	// Add ixn to lattice
	latt->all_conns_vh_set_ixn_dict(ixn_dict);

	// Add all units to moment of ixn func
	latt->all_units_h_add_to_moment_b(ixn_bX->get_moment());
	latt->all_conns_vh_add_to_moment_w(ixn_wAX->get_moment());

	cout << " > end hidden" << endl;

	cout << "--- [Finished] Making lattice ---" << endl;
	cout << endl;

	/****************************************
	Do an iteration of the learning problem
	****************************************/

	// Opt solver class
	OptProblem opt(latt);

	// Add ixn params
	opt.add_ixn_param(ixn_hA);
	opt.add_ixn_param(ixn_bX);
	opt.add_ixn_param(ixn_wAX);

	// Params
	int no_opt_steps = 100;
	int no_timesteps = 100;
	double dt=0.1;
	int no_latt_sampling_steps = 10;
	int batch_size = 5;
	double dopt=0.001;

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

		// Adjust integral range
		/*
		if (opt_step <= 50) {
			options.VAL_random_integral_range_size = 10;
		} else if (opt_step <= 75) {
			options.VAL_random_integral_range_size = 25;
		} else {
			options.VAL_random_integral_range_size = 50;
		};
		*/

		// Take a step
		opt.solve_one_step(opt_step,no_timesteps,batch_size,dt,dopt,no_latt_sampling_steps,fnames,options);

		// Sample moments (whole trajectory)
		if (opt_step==no_opt_steps) {
			opt.wake_sleep_loop(0,no_timesteps,batch_size,no_latt_sampling_steps,fnames,false);
		};

		// Print moments
		std::cout << ixn_hA->get_name() << " " << std::flush;
		ixn_hA->get_moment()->print_moment_comparison();

		std::cout << ixn_wAX->get_name() << " " << std::flush;
		ixn_wAX->get_moment()->print_moment_comparison();

		std::cout << ixn_bX->get_name() << " " << std::flush;
		ixn_bX->get_moment()->print_moment_comparison();

		// Write moments
		ixn_hA->get_moment()->write_to_file("data_learned/moments/"+ixn_hA->get_name()+"_"+pad_str(opt_step,3)+".txt");
		ixn_wAX->get_moment()->write_to_file("data_learned/moments/"+ixn_wAX->get_name()+"_"+pad_str(opt_step,3)+".txt");
		ixn_bX->get_moment()->write_to_file("data_learned/moments/"+ixn_bX->get_name()+"_"+pad_str(opt_step,3)+".txt");

		// Write ixn params
		ixn_hA->write_to_file("data_learned/ixn_params/"+ixn_hA->get_name()+"_"+pad_str(opt_step,3)+".txt");
		ixn_wAX->write_to_file("data_learned/ixn_params/"+ixn_wAX->get_name()+"_"+pad_str(opt_step,3)+".txt");
		ixn_bX->write_to_file("data_learned/ixn_params/"+ixn_bX->get_name()+"_"+pad_str(opt_step,3)+".txt");

		// Write adjoint
		ixn_hA->get_adjoint()->write_to_file("data_learned/adjoint/"+ixn_hA->get_name()+"_"+pad_str(opt_step,3)+".txt");
		ixn_wAX->get_adjoint()->write_to_file("data_learned/adjoint/"+ixn_wAX->get_name()+"_"+pad_str(opt_step,3)+".txt");
		ixn_bX->get_adjoint()->write_to_file("data_learned/adjoint/"+ixn_bX->get_name()+"_"+pad_str(opt_step,3)+".txt");

		// Write diff eq rhs
		if (opt_step == no_opt_steps) {
			rhs_hA->write_to_file("data_learned/diff_eq_rhs/"+ixn_hA->get_name()+"_"+pad_str(opt_step,3)+".txt");
			rhs_wAX->write_to_file("data_learned/diff_eq_rhs/"+ixn_wAX->get_name()+"_"+pad_str(opt_step,3)+".txt");
			rhs_bX->write_to_file("data_learned/diff_eq_rhs/"+ixn_bX->get_name()+"_"+pad_str(opt_step,3)+".txt");
		};
	};

	return 0;
};