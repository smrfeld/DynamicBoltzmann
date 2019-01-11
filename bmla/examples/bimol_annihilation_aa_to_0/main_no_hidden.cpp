#include <bmla>

#include <iostream>

using namespace bmla;
using namespace std;

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

	auto ixn_hA = make_shared<IxnParam>("hA",IxnParamType::H, 0.0);

	cout << "--- [Finished] Making ixn func ---" << endl;
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

	/****************************************
	Do an iteration of the learning problem
	****************************************/

	// Opt solver class
	OptProblem opt(latt);

	// Add ixn params
	opt.add_ixn_param(ixn_hA);

	// Params
	int no_opt_steps = 1000;
	int no_latt_sampling_steps = 10;
	int batch_size = 5;
	double dopt=0.0001;

	// Filenames
	FNameColl fnames;
	for (auto i_batch=0; i_batch<10; i_batch++) {
		// Add
		fnames.add_fname("stoch_sim/lattice_v"+pad_str(i_batch,3)+"/lattice/0000.txt");
	};

	// Options
	OptionsSolve options;
	options.VERBOSE_WAKE_ASLEEP = false;
	options.VERBOSE_MOMENT = false;

	// Write ixn params
	ixn_hA->write_to_file("data_learned_no_hidden/ixn_params/"+ixn_hA->get_name()+"_"+pad_str(0,3)+".txt");

	// Manually
	opt.init_structures(batch_size);
	for (auto opt_step=1; opt_step<=no_opt_steps; opt_step++) {

		std::cout << "------------------" << std::endl;
		std::cout << "Opt step: " << opt_step << " / " << no_opt_steps << std::endl;
		std::cout << "------------------" << std::endl;

		// Take a step
		opt.solve_one_step(batch_size,dopt,no_latt_sampling_steps,fnames,options);

		// Sample moments (whole trajectory)
		if (opt_step==no_opt_steps) {
			opt.wake_sleep_loop(batch_size,no_latt_sampling_steps,fnames,false);
		};

		// Print moments
		std::cout << ixn_hA->get_name() << " " << ixn_hA->get_val() << std::endl;
		// ixn_hA->get_moment()->print_moment_comparison();

		// Write moments
		ixn_hA->get_moment()->write_to_file("data_learned_no_hidden/moments/"+ixn_hA->get_name()+"_"+pad_str(opt_step,3)+".txt");

		// Write ixn params
		ixn_hA->write_to_file("data_learned_no_hidden/ixn_params/"+ixn_hA->get_name()+"_"+pad_str(opt_step,3)+".txt");
	};

	return 0;
};