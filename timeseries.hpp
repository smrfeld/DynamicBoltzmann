/************************************
* Namespace for DynamicBoltzmann
************************************/

namespace DynamicBoltzmann {

	/****************************************
	Timeseries
	****************************************/

	class Timeseries
	{
	private:

		// Range in time
		double _tmin;
		double _tmax;

		// Number timesteps
		int _n_t;

		// Timestep
		double _dt;

		// Grid to store solution
		// Pointer so size can be dynamically allocated by constructor
		double* _soln;

	public:

		// Constructor
		Timeseries(double tmin, double tmax, int n_t);
		Timeseries(const Timeseries& t);
		Timeseries& operator=(const Timeseries& t);
		~Timeseries();

	};

};