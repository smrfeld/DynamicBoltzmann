/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	HiddenLayer
	****************************************/

	/********************
	Constructor
	********************/

	HiddenLayer::HiddenLayer() {};
	HiddenLayer::HiddenLayer(const HiddenLayer& other) {
		_copy(other);
	};
	HiddenLayer::HiddenLayer(HiddenLayer&& other) {
		_copy(other);
		other._reset();
	};
	HiddenLayer& HiddenLayer::operator=(const HiddenLayer& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	HiddenLayer& HiddenLayer::operator=(HiddenLayer&& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
	};
	HiddenLayer::~HiddenLayer() {
		_clean_up();
	};
	void HiddenLayer::_clean_up() {
		// Nothing....
	};
	void HiddenLayer::_reset() {
		_layer.clear();
		_sp_map.clear();
		_sp_vec.clear();
	};
	void HiddenLayer::_copy(const HiddenLayer& other) {
		_layer = other._layer;
		_sp_map = other._sp_map;
		_sp_vec = other._sp_vec;	
	};

	/********************
	Add a species
	********************/

	void HiddenLayer::add_species_possibility(Species *sp) {
		_sp_map[sp->name()] = sp;
		_sp_vec.push_back(sp);
	};

	/********************
	Activate hidden layer
	********************/

	void activate(bool binary=true) {

		if (!binary) {
			std::cerr << "ERROR: Not binary currently unsuportted" << std::endl;
			exit(EXIT_FAILURE);
		};

		// Propensities
		std::vector<double> props;

		// Go through all units
		for (hidden_layer_it hit = _layer.begin(); hit != _layer.end(); hit++) {

			// Clear props, probs
			props.clear();
			props.push_back(0.0);

			// Empty = 1
			props.push_back(1.0);

			// Go through all possible species this could be, calculate propensities
			for (auto sp_new: _sp_vec) {

				// Bias
				energy = sp_new->b();

				// Go through connections to this layer
				energy += hit->get_activation(sp_new);

				// Append prop
				props.push_back(props.back()+exp(energy));
			};
		};

		// Commit probs
		if (binary) {

			// Sample RV
			i_chosen = sample_prop_vec(props);

			if (i_chosen==0) {
				// Flip down (new spin = 0)
				it->set_site_empty();
			} else {
				// Make the appropriate species at this site (guaranteed empty)
				it->set_site_binary(_sp_vec[i_chosen-1]);
			};
		};
	};

};