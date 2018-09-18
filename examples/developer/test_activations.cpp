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

	double init_val = 0.8;
	auto ixn = make_shared<IxnParam>("bias for A",IxnParamType::H, init_val);
	ixn->add_species_h(species_A);

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
	Get moment
	****************************************/

	auto moment = ixn->get_moment();

	cout << "Got moment of " << ixn->get_name() << endl;
	cout << endl;

	/****************************************
	Make visible unit
	****************************************/

	// Possible species
	vector<Sptr> species_possible({species_A});

	// Units
	UnitVisible vis_unit_1 = UnitVisible(1,species_possible);

	// Make the occupied species of the units each be A
	vis_unit_1.set_binary(species_A);

	cout << "Made visible unit" << endl;
	cout << endl;

	/****************************************
	Make unit ixn (bias) dict
	****************************************/

	auto bias_dict = make_shared<BiasDict>();
	bias_dict->add_ixn(species_A,ixn);

	// Add to unit 1
	vis_unit_1.set_bias_dict(bias_dict);

	cout << "Added bias to unit 1 for species A" << endl;
	cout << endl;

	/****************************************
	Make another visible unit
	****************************************/

	auto species_B = make_shared<Species>("B");
	species_possible = {species_A, species_B};
	UnitVisible vis_unit_2 = UnitVisible(2, species_possible);
	UnitVisible vis_unit_3 = UnitVisible(3, species_possible);

	// Occupy third unit
	vis_unit_3.set_binary(species_A);

	// Add bias to unit 2
	vis_unit_2.set_bias_dict(bias_dict);

	cout << "Made 2 more visible units" << endl;
	cout << endl;

	/****************************************
	Make some order 2 ixns
	****************************************/

	// A-A
	auto ixn_AA = make_shared<IxnParam>("A-A ixn",IxnParamType::J, init_val);
	ixn_AA->add_species_j(species_A,species_A);

	// A-B
	auto ixn_AB = make_shared<IxnParam>("A-B ixn",IxnParamType::J, 0.3);
	ixn_AB->add_species_j(species_A,species_B);

	cout << "Made order 2 ixns: " << ixn_AA->get_name() << " and " << ixn_AB->get_name() << endl;
	cout << endl;

	/****************************************
	Make a connection, add to units
	****************************************/

	// Connections
	ConnVV conn_12 = ConnVV(&vis_unit_1,&vis_unit_2);
	vis_unit_1.add_conn(&conn_12,0);
	vis_unit_2.add_conn(&conn_12,1);
	ConnVV conn_23 = ConnVV(&vis_unit_2,&vis_unit_3);
	vis_unit_2.add_conn(&conn_23,0);
	vis_unit_3.add_conn(&conn_23,1);

	// Ixns that can live on these connection
	auto o2_dict = make_shared<O2IxnDict>();
	o2_dict->add_ixn(species_A,species_A,ixn_AA);
	o2_dict->add_ixn(species_A,species_B,ixn_AB);
	/*
	o2_dict->add_ixn_with_all_permutations(species_A,species_A,ixn_AA);
	o2_dict->add_ixn_with_all_permutations(species_A,species_B,ixn_AB);
	*/

	// Add to the conns
	conn_12.set_ixn_dict(o2_dict);
	conn_23.set_ixn_dict(o2_dict);

	cout << "Made connections from unit 1->2->3" << endl;
	cout << endl;

	/****************************************
	Get activation/sample
	****************************************/

	// Species A @ unit 2 

	double act_species_A_at_unit_2 = vis_unit_2.get_act_for_species_at_timepoint(species_A,0);
	cout << "Activation of species A at unit 2: " << act_species_A_at_unit_2 << endl;

	double act_species_A_at_unit_2_from_bias = bias_dict->get_ixn_at_timepoint(species_A,0);
	cout << "    from bias: " << act_species_A_at_unit_2_from_bias << endl;

	double act_species_A_at_unit_2_from_conn_12 = conn_12.get_act_for_species_at_unit_at_timepoint(species_A,1,0);
	cout << "    from 1->2 connection: " << act_species_A_at_unit_2_from_conn_12 << endl;

	double act_species_A_at_unit_2_from_conn_23 = conn_23.get_act_for_species_at_unit_at_timepoint(species_A,0,0);
	cout << "    from 2->3 connection: " << act_species_A_at_unit_2_from_conn_23 << endl;

	// Species B @ unit 2

	double act_species_B_at_unit_2 = vis_unit_2.get_act_for_species_at_timepoint(species_B,0);
	cout << "Activation of species B at unit 2: " << act_species_B_at_unit_2 << endl;

	double act_species_B_at_unit_2_from_bias = bias_dict->get_ixn_at_timepoint(species_B,0);
	cout << "    from bias: " << act_species_B_at_unit_2_from_bias << endl;

	double act_species_B_at_unit_2_from_conn_12 = conn_12.get_act_for_species_at_unit_at_timepoint(species_B,1,0);
	cout << "    from 1->2 connection: " << act_species_B_at_unit_2_from_conn_12 << endl;

	double act_species_B_at_unit_2_from_conn_23 = conn_23.get_act_for_species_at_unit_at_timepoint(species_B,0,0);
	cout << "    from 2->3 connection: " << act_species_B_at_unit_2_from_conn_23 << endl;

	// Sample unit 2

	double prob_0, prob_A, prob_B;
	for (auto i=0; i<10; i++)
	{
		vis_unit_2.sample_at_timepoint(0);
		prob_0 = vis_unit_2.get_prob(nullptr);
		prob_A = vis_unit_2.get_prob(species_A);
		prob_B = vis_unit_2.get_prob(species_B);
		cout << "Sampled unit 2: now prob empty = " << prob_0 << " A = " << prob_A << " B = " << prob_B << endl;
	};

	return 0;
};




