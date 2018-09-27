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

	auto ixn_hA = make_shared<IxnParam>("hA",IxnParamType::H, -0.1);
	auto ixn_wAX = make_shared<IxnParam>("wAX",IxnParamType::W, 0.1);
	auto ixn_bX = make_shared<IxnParam>("bX",IxnParamType::B, -0.1);

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
	TEMP: TEST HIDDEN LAYER
	****************************************/

	/*
	// Convert to binary
	latt->all_units_convert_to_b_mode();

	// Read in some file
	latt->read_from_file("stoch_sim/lattice_v000/lattice/0000.txt");

	// Activate hidden
	latt->sample_h();

	// Print hiddens
	int ctr_x=0;
	for (auto i=1; i<1000; i++) { // 999 total
		Sptr sp = latt->get_unit_h(layer,i).get_b_mode_species();
		if (sp) {
			ctr_x++;
			// cout << i << " " << sp->get_name() << endl;
		} else {
			// cout << i << " " << "empty" << endl;
		};
	};
	cout << "X: " << ctr_x << endl;

	// Convert to prob
	latt->all_units_convert_to_p_mode();

	// Activate visible
	latt->sample_v(false);
	latt->sample_h(false);
	for (auto i=1; i<=1000; i++) {
		latt->get_unit_v(i).print();
	};
	for (auto i=1; i<1000; i++) {
		latt->get_unit_h(layer,i).print();
	};

	return 0;
	*/

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
	int no_opt_steps = 1000;
	int no_latt_sampling_steps = 10;
	int batch_size = 5;
	double dopt=0.001;

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
	ixn_hA->write_to_file("data_learned/ixn_params/"+ixn_hA->get_name()+"_"+pad_str(0,3)+".txt");
	ixn_wAX->write_to_file("data_learned/ixn_params/"+ixn_wAX->get_name()+"_"+pad_str(0,3)+".txt");
	ixn_bX->write_to_file("data_learned/ixn_params/"+ixn_bX->get_name()+"_"+pad_str(0,3)+".txt");

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

		std::cout << ixn_wAX->get_name() << " " << ixn_wAX->get_val() << std::endl;
		// ixn_wAX->get_moment()->print_moment_comparison();

		std::cout << ixn_bX->get_name() << " " << ixn_bX->get_val() << std::endl;
		// ixn_bX->get_moment()->print_moment_comparison();

		// Write moments
		ixn_hA->get_moment()->write_to_file("data_learned/moments/"+ixn_hA->get_name()+"_"+pad_str(opt_step,3)+".txt");
		ixn_wAX->get_moment()->write_to_file("data_learned/moments/"+ixn_wAX->get_name()+"_"+pad_str(opt_step,3)+".txt");
		ixn_bX->get_moment()->write_to_file("data_learned/moments/"+ixn_bX->get_name()+"_"+pad_str(opt_step,3)+".txt");

		// Write ixn params
		ixn_hA->write_to_file("data_learned/ixn_params/"+ixn_hA->get_name()+"_"+pad_str(opt_step,3)+".txt");
		ixn_wAX->write_to_file("data_learned/ixn_params/"+ixn_wAX->get_name()+"_"+pad_str(opt_step,3)+".txt");
		ixn_bX->write_to_file("data_learned/ixn_params/"+ixn_bX->get_name()+"_"+pad_str(opt_step,3)+".txt");
	};

	return 0;
};