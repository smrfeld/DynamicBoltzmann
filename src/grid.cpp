#include "grid.hpp"

// Other headers
#include "../include/dynamicboltz_bits/general.hpp"

#include <iostream>
#include "math.h"
#include <fstream>

/************************************
* Namespace for dboltz
************************************/

namespace dboltz {

	/****************************************
	Grid
	****************************************/

	/********************
	Constructor
	********************/

	Grid::Grid(std::string name, double min, double max, int n)
	{
		_name = name;
		_min = min;
		_max = max;
		_n = n;
		_delta = (max-min)/(n-1);
		_grid = new double[n];
		for (int i=0; i<n; i++) {
			_grid[i] = min+i*_delta;
		};
	};
	Grid::Grid(const Grid& other) 
	{
		_copy(other);
	};
	Grid& Grid::operator=(const Grid& other)
	{
		if (this != &other)
		{			
			_clean_up();
			_copy(other);
		};
		return *this;		
	};
	Grid::~Grid() {
		_clean_up();
	};
	void Grid::_copy(const Grid& other) {
		_name = other._name;
		_min = other._min;
		_max = other._max;
		_n = other._n;
		_delta = other._delta;
		_grid = new double[_n];
		std::copy( other._grid, other._grid + _n, _grid );
	};
	void Grid::_clean_up() {
		safeDelArr(_grid);
	};
	
	/********************
	Getters/setters
	********************/

	std::string Grid::name() const { return _name; };

	double Grid::delta() const { return _delta; };

	int Grid::n() const { return _n; };

	double Grid::get_by_idx(int i) const {
		return _grid[i];
	};

	/********************
	Set n
	********************/

	// Set n
	void Grid::set_n(int n) {
		// Clean old
		safeDelArr(_grid);
		// Set new
		_n = n;
		// New grid
		_grid = new double[_n];
		for (int i=0; i<_n; i++) {
			_grid[i] = _min+i*_delta;
		};
	};

	/********************
	Surrounding idxs
	********************/

	int Grid::surrounding_idxs(double x) const
	{
		int i = (x - _min) / _delta;
		if (i==_n-1) {i--;};
		return i;
	};

	double Grid::frac_between(double x) const
	{
		return frac_between(x,surrounding_idxs(x));
	};

	double Grid::frac_between(double x, int i) const
	{
		return (x-_grid[i]) / _delta;
	};

	/********************
	Check if a given point is in the grid
	********************/

	bool Grid::in_grid(double x) const
	{
		if (x < _min || x > _max) { 
			return false; 
		} else {
			return true;
		};
	};
 
	/********************
	Print grid range
	********************/

	void Grid::print_grid_range() const
	{
		std::cout << "Grid: " << _name << " min: " << _min << " max: " << _max << std::endl;
	};

	/********************
	Write grid into an ofstream
	********************/

	void Grid::write_grid(std::string fname) const {
		std::ofstream f;
		f.open(fname);
		for (int i=0; i<_n; i++) {
			f << _grid[i] << "\n";
		};
		f.close();
	};

	/********************
	Test: if a sin func were defined on the grid
	********************/

	std::vector<double> Grid::test_sin() const {
		std::vector<double> x;
		for (int i=0; i<_n; i++) {
			x.push_back(sin(2*M_PI*_grid[i]/(_max-_min)));
		};
		return x;
	};
	std::vector<double> Grid::test_cos() const {
		std::vector<double> x;
		for (int i=0; i<_n; i++) {
			x.push_back(cos(2*M_PI*_grid[i]/(_max-_min)));
		};
		return x;	
	};

};