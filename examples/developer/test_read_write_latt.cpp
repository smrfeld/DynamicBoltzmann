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
	Make lattice
	****************************************/

	Lattice latt(1,10);

	// Set possible species
	latt.all_units_v_add_possible_species(species_A);
	latt.all_units_v_add_possible_species(species_B);

	/****************************************
	Populate the lattice
	****************************************/

	// Populate randomly
	latt.get_unit_v(3).set_b_mode_species(species_A);
	latt.get_unit_v(5).set_b_mode_species(species_B);
	latt.get_unit_v(7).set_b_mode_species(species_A);

	cout << "Made some sites occupied" << endl;
	latt.print_occupancy();
	cout << endl;

	/****************************************
	Write
	****************************************/

	latt.write_to_file("test_read_write_latt_1.txt");

	cout << "Wrote latt" << endl;
	cout << endl;

	/****************************************
	Read
	****************************************/

	latt.read_from_file("test_read_write_latt_2.txt");

	cout << "Read different latt" << endl;
	latt.print_occupancy();
	cout << endl;



	return 0;
};




