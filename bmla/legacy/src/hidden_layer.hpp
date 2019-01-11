/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	typedef std::list<HiddenUnit> hidden_layer;
	typedef std::list<HiddenUnit>::iterator hidden_layer_it;

	/****************************************
	HiddenLayer
	****************************************/

	class HiddenLayer
	{
	private:

		// Internal maps
		hidden_layer _layer;

		// Pointers to species present
		std::map<std::string,HiddenSpecies*> _sp_map;
		std::vector<HiddenSpecies*> _sp_vec;

		// Contructor helpers
		void _clean_up();
		void _reset();
		void _copy(const Lattice& other);

	public:

		/********************
		Constructor
		********************/

		HiddenLayer();
		HiddenLayer(const Lattice& other);
		HiddenLayer(HiddenLayer&& other);
		HiddenLayer& operator=(const HiddenLayer& other);
		HiddenLayer& operator=(HiddenLayer&& other);
		~HiddenLayer();

		/********************
		Add a species
		********************/

		void add_species_possibility(HiddenSpecies *sp);

		/********************
		Activate hidden layer
		********************/

		void activate(bool binary=true);

	};

};