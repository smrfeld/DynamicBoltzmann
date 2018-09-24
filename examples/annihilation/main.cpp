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
	Lattice latt(3,10);

	// Set possible species
	latt.all_units_v_add_possible_species(species_A);

	// Make ixn dict
	auto bias_dict = make_shared<BiasDict>();
	bias_dict->add_ixn(species_A,ixn);

	// Add ixn to lattice
	latt.all_units_v_set_bias_dict(bias_dict);

	// Add all units to moment of ixn func
	for (auto i1=1; i1<= 10; i1++) {
		for (auto i2=1; i2<= 10; i2++) {
			for (auto i3=1; i3<= 10; i3++) {
				ixn->get_moment()->add_unit_to_monitor_h(&latt.get_unit_v(i1,i2,i3));
			};
		};
	};

	cout << "--- [Finished] Making lattice ---" << endl;
	cout << endl;

	/****************************************
	Do an iteration of the learning problem
	****************************************/

	/********************
	Params
	********************/

	int no_opt_steps = 50;
	int no_timesteps = 20;
	double dt=0.1;
	int no_latt_sampling_steps = 10;
	int batch_size = 5;
	double dopt=0.002;

	/********************
	Verbose options
	********************/

	bool VERBOSE_NU = false;
	bool VERBOSE_ADJOINT = false;
	bool VERBOSE_UPDATE = false;
	bool VERBOSE_SAMPLE = false;
	bool VERBOSE_MOMENT = true;

	/********************
	Init structures
	********************/

	// Ixn func
	ixn->set_no_timesteps(no_timesteps);

	// Adjoint
	adjoint->set_no_timesteps(no_timesteps);
	adjoint->set_zero_end_cond_timepoint(no_timesteps);

	// Moment
	ixn->get_moment()->set_no_timesteps(no_timesteps);
	ixn->get_moment()->set_batch_size(batch_size);

	/********************
	Go...
	********************/

	for (int i_opt_step=1; i_opt_step<=no_opt_steps; i_opt_step++)
	{

		std::cout << "------------------" << std::endl;
		std::cout << "Opt step: " << i_opt_step << " / " << no_opt_steps << std::endl;
		std::cout << "------------------" << std::endl;

		/********************
		Solve diff eq for F
		********************/

		if (VERBOSE_NU) {
			cout << "--- Solving diff eq ---" << endl;
			cout << "	t = " << 0 << " val = " << ixn->get_val_at_timepoint(0) << endl;
		};
		for (auto t=0; t<no_timesteps; t++) {
			ixn->solve_diff_eq_at_timepoint_to_plus_one(t,dt);

			if (VERBOSE_NU) {
				cout << "	t = " << t+1 << " val = " << ixn->get_val_at_timepoint(t+1) << endl;
			};
		};

		if (VERBOSE_NU) {
			cout << "--- [Finished] Solving diff eq ---" << endl;
			cout << endl;
		};

		/********************
		Sample lattice
		********************/

		if (VERBOSE_SAMPLE) {
			cout << "--- Sampling lattice ---" << endl;
		};

		for (int timepoint=0; timepoint<=no_timesteps; timepoint++) {

			if (VERBOSE_SAMPLE) {
				cout << "(" << timepoint << "/" << no_timesteps << ") : " << flush;
			};

			for (int i_batch=0; i_batch<batch_size; i_batch++)
			{

				if (VERBOSE_SAMPLE) {
					cout << "." << flush;
				};

				// Read latt
				latt.read_from_file("stoch_sim/lattice_v"+pad_str(randI(0,9),2)+"/lattice/"+pad_str(timepoint,4)+".txt");

				// Reap awake
				ixn->get_moment()->reap_as_timepoint_in_batch(MomentType::AWAKE, timepoint, i_batch);

				// Sample
				for (int i_sampling_step=0; i_sampling_step<no_latt_sampling_steps; i_sampling_step++) 
				{
					latt.sample_v_at_timepoint(timepoint);
				};

				// Print
				// latt.print_occupancy();

				// Reap asleep
				ixn->get_moment()->reap_as_timepoint_in_batch(MomentType::ASLEEP, timepoint, i_batch);
			};
			
			// Average moments at this timepoint
			ixn->get_moment()->average_reaps_as_timepoint(MomentType::AWAKE, timepoint);
			ixn->get_moment()->average_reaps_as_timepoint(MomentType::ASLEEP, timepoint);

			if (VERBOSE_SAMPLE) {
				cout << endl;
			};

		};

		if (VERBOSE_MOMENT) {
			ixn->get_moment()->print_moment_comparison();
		};

		if (VERBOSE_SAMPLE) {
			cout << "--- [Finished] Sampled lattice ---" << endl;
			cout << endl;
		};

		/********************
		Solve diff eq for adjoint
		********************/

		if (VERBOSE_ADJOINT) {
			cout << "--- Solving adjoint ---" << endl;
			cout << "	t = " << no_timesteps << " val = " << adjoint->get_val_at_timepoint(no_timesteps) << endl;
		};
		for (auto t=no_timesteps; t>=1; t--) {
			adjoint->solve_diff_eq_at_timepoint_to_minus_one(t,dt);

			if (VERBOSE_ADJOINT) {
				cout << "	t = " << t-1 << " val = " << adjoint->get_val_at_timepoint(t-1) << endl;
			};
		};

		if (VERBOSE_ADJOINT) {
			cout << "--- [Finished] Solving adjoint ---" << endl;
			cout << endl;
		};

		/********************
		Form the update
		********************/

		if (VERBOSE_UPDATE) {
			cout << "--- Calculating update ---" << endl;
		};

		// rhs->calculate_update(0,no_timesteps,dt,dopt);
		rhs->calculate_update(0,randI(5,20),dt,dopt);

		if (VERBOSE_UPDATE) {
			rhs->print_update_stored();

			cout << "--- [Finished] Calculating update ---" << endl;
			cout << endl;
		};

		/********************
		Committ the update
		********************/

		if (VERBOSE_UPDATE) {
			cout << "--- Committing update ---" << endl;
		};

		rhs->committ_update();

		if (VERBOSE_UPDATE) {
			cout << "New RHS grid:" << endl;
			rhs->print_grid_pts_inside();

			cout << "--- [Finished] Committing update ---" << endl;
			cout << endl;
		};

	};

	return 0;
};