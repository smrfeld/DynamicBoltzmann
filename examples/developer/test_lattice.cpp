#include <dynamicboltz>

#include <iostream>

using namespace dblz;
using namespace std;

int main() {

	/****************************************
	Make some species
	****************************************/

	auto species_A = make_shared<Species>("A");
	auto species_B = make_shared<Species>("B");

	cout << "Made species: " << species_A->get_name() << endl;
	cout << "Made species: " << species_B->get_name() << endl;
	cout << endl;

	/****************************************
	Make ixn
	****************************************/

	auto ixn_bias_A = make_shared<IxnParam>("bias for A",IxnParamType::H, 0.8);
	ixn_bias_A->add_species_h(species_A);

	// A-A
	auto ixn_AA = make_shared<IxnParam>("A-A ixn",IxnParamType::J, 0.8);
	ixn_AA->add_species_j(species_A,species_A);

	// A-B
	auto ixn_AB = make_shared<IxnParam>("A-B ixn",IxnParamType::J, 0.3);
	ixn_AB->add_species_j(species_A,species_B);

	cout << "Made ixn bias for A: " << ixn_bias_A->get_name() << endl;
	cout << "Made ixn A-A: " << ixn_AA->get_name() << endl;
	cout << "Made ixn A-B: " << ixn_AB->get_name() << endl;
	cout << endl;

	/****************************************
	Make lattice
	****************************************/

	Lattice latt(1,10);

	// Set possible species
	latt.add_possible_species_to_all_units_visible(species_A);
	latt.add_possible_species_to_all_units_visible(species_B);

	// Make NN connectivity
	latt.init_conns_NN_all_units_visible();

	cout << "Made lattice" << endl;
	cout << endl;

	/****************************************
	Set ixn dicts
	****************************************/

	auto ixn_dict = make_shared<O2IxnDict>();
	ixn_dict->add_ixn_and_its_reverse(species_A,species_A,ixn_AA);
	ixn_dict->add_ixn_and_its_reverse(species_A,species_B,ixn_AB);

	auto bias_dict = make_shared<BiasDict>();
	bias_dict->add_ixn(species_A,ixn_bias_A);

	// Add to lattice
	latt.set_bias_dict_of_all_units_visible(bias_dict);
	latt.set_ixn_dict_all_conns_vv(ixn_dict);

	cout << "Added ixns to lattice conns" << endl;
	cout << endl;

	/****************************************
	Populate the lattice
	****************************************/

	// Populate randomly
	latt.get_unit_visible(3).set_b_mode_species(species_A);
	latt.get_unit_visible(5).set_b_mode_species(species_B);
	latt.get_unit_visible(7).set_b_mode_species(species_A);

	cout << "Made some sites occupied" << endl;
	latt.print_occupancy();
	cout << endl;

	/****************************************
	Sample
	****************************************/

	for (auto i=0; i<10; i++) {
		latt.sample_at_timepoint(0);

		cout << "Sampled:" << endl;
		cout << "Count A = " << latt.get_count(species_A) << endl;
		cout << "Count B = " << latt.get_count(species_B) << endl;
		cout << "NN(AA) = " << latt.get_count(species_A,species_A) << endl;
		cout << "NN(AB - this dir) = " << latt.get_count(species_A,species_B,true) << endl;
		cout << "NN(BA - this dir) = " << latt.get_count(species_A,species_B,true) << endl;
		cout << "NN(AB - both dir) = " << latt.get_count(species_A,species_B) << endl;
		cout << "NN(BB) = " << latt.get_count(species_B,species_B) << endl;
		cout << "Lattice:" << endl;
		latt.print_occupancy();
		cout << endl;
	};

	/****************************************
	Add units to the ixn funcs moment
	****************************************/

	auto moment_A = ixn_bias_A->get_moment();
	moment_A->set_batch_size(10);

	auto moment_AB = ixn_AB->get_moment();
	moment_AB->set_batch_size(10);

	// Add
	for (auto i=1; i<=10; i++) {
		moment_A->add_unit_to_monitor_h(&latt.get_unit_visible(i));
	};
	for (auto i=1; i<=9; i++) {
		moment_AB->add_conn_to_monitor_j(&latt.get_conn_vv(i,i+1));
	};

	cout << "Setup moment for ixn func bias A and A-B" << endl;
	cout << endl;

	/****************************************
	Reap moment
	****************************************/

	for (auto i=0; i<10; i++) {
		latt.sample_at_timepoint(0);
		cout << "Count A = " << latt.get_count(species_A) << " NN(AB) = " << latt.get_count(species_A,species_B) << endl;
		moment_A->reap_as_timepoint_in_batch(MomentType::AWAKE, 0, i);
		moment_AB->reap_as_timepoint_in_batch(MomentType::AWAKE, 0, i);
	};
	moment_A->average_reaps_as_timepoint(MomentType::AWAKE, 0);
	moment_AB->average_reaps_as_timepoint(MomentType::AWAKE, 0);

	cout << "Reaped and averaged moment while sampling: bias A = " << moment_A->get_moment_at_timepoint(MomentType::AWAKE,0) << " NN(AB) = " << moment_AB->get_moment_at_timepoint(MomentType::AWAKE,0) << endl;
	cout << endl;

	return 0;
};




