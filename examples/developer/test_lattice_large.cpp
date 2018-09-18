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

	// A-A
	auto ixn_AA = make_shared<IxnParam>("A-A ixn",IxnParamType::J, 0.8);

	// A-B
	auto ixn_AB = make_shared<IxnParam>("A-B ixn",IxnParamType::J, 0.3);

	cout << "Made ixn bias for A: " << ixn_bias_A->get_name() << endl;
	cout << "Made ixn A-A: " << ixn_AA->get_name() << endl;
	cout << "Made ixn A-B: " << ixn_AB->get_name() << endl;
	cout << endl;

	/****************************************
	Make lattice
	****************************************/
	
	// 20 x 20 x 20
	Lattice latt(3,20);

	// Set possible species
	latt.add_possible_species_to_all_units_vis(species_A);
	latt.add_possible_species_to_all_units_vis(species_B);

	// Make NN connectivity
	latt.init_conns_NN_all_units_vis();

	cout << "Made lattice" << endl;
	cout << endl;

	/****************************************
	Set ixn dicts
	****************************************/

	auto ixn_dict = make_shared<O2IxnDict>();
	ixn_dict->add_ixn(species_A,species_A,ixn_AA);
	ixn_dict->add_ixn(species_A,species_B,ixn_AB,false); // one-way, not reversibly

	auto bias_dict = make_shared<BiasDict>();
	bias_dict->add_ixn(species_A,ixn_bias_A);

	// Add to lattice
	latt.set_bias_dict_of_all_units_vis(bias_dict);
	latt.set_ixn_dict_of_all_conns_vv(ixn_dict);

	cout << "Added ixns to lattice conns" << endl;
	cout << endl;

	/****************************************
	Add units to the ixn funcs moment
	****************************************/

	auto moment_A = ixn_bias_A->get_moment();
	moment_A->set_batch_size(10);

	auto moment_AB = ixn_AB->get_moment();
	moment_AB->set_batch_size(10);

	auto moment_AA = ixn_AA->get_moment();
	moment_AA->set_batch_size(10);

	// Add
	latt.add_all_units_vis_to_moment_h(moment_A);
	latt.add_all_conns_vv_to_moment_j(moment_AA);
	latt.add_all_conns_vv_to_moment_j(moment_AB);

	cout << "Setup moment for ixn func bias A and A-B and A-A" << endl;
	cout << endl;

	/****************************************
	Reap moment
	****************************************/

	for (auto i=0; i<10; i++) {
		latt.sample_at_timepoint(0);
		// cout << "Count A = " << latt.get_count(species_A) << " NN(AA) = " << latt.get_count(species_A,species_A) << " NN(AB) = " << latt.get_count(species_A,species_B) << endl;
		moment_A->reap_as_timepoint_in_batch(MomentType::AWAKE, 0, i);
		moment_AB->reap_as_timepoint_in_batch(MomentType::AWAKE, 0, i);
		moment_AA->reap_as_timepoint_in_batch(MomentType::AWAKE, 0, i);
	};
	moment_A->average_reaps_as_timepoint(MomentType::AWAKE, 0);
	moment_AB->average_reaps_as_timepoint(MomentType::AWAKE, 0);
	moment_AA->average_reaps_as_timepoint(MomentType::AWAKE, 0);

	cout << "Reaped and averaged moment while sampling: bias A = " << moment_A->get_moment_at_timepoint(MomentType::AWAKE,0) << " NN(AA) = " << moment_AA->get_moment_at_timepoint(MomentType::AWAKE,0) << " NN(AB) = " << moment_AB->get_moment_at_timepoint(MomentType::AWAKE,0) << endl;
	cout << endl;

	return 0;
};




