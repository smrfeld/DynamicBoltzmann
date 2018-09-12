#include "domain.hpp"

// Other headers
#include "../include/dynamicboltz_bits/general.hpp"

#include <iostream>
#include "math.h"
#include <fstream>

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	Grid
	****************************************/

	/********************
	Constructor
	********************/

	Domain1D::Domain1D(std::string name, double min, double max, int no_pts)
	{
		_name = name;
		_min = min;
		_max = max;
		_no_pts = no_pts;
		_delta = (_max-_min)/(_no_pts-1);
	};
	Domain1D::Domain1D(const Domain1D& other) 
	{
		_copy(other);
	};
	Domain1D& Domain1D::operator=(const Domain1D& other)
	{
		if (this != &other)
		{			
			_clean_up();
			_copy(other);
		};
		return *this;		
	};
	Domain1D::~Domain1D() {
		_clean_up();
	};
	void Domain1D::_copy(const Domain1D& other) {
		_name = other._name;
		_min = other._min;
		_max = other._max;
		_no_pts = other._no_pts;
		_delta = other._delta;
	};
	void Domain1D::_clean_up() {
	};
	
	/********************
	Getters
	********************/

	std::string Domain1D::get_name() const { return _name; };

	double Domain1D::get_delta() const { return _delta; };

	int Domain1D::get_no_pts() const { return _no_pts; };

	/********************
	Value of a point
	********************/

	double Domain1D::get_by_idx(int i) const {
		if (i >= _no_pts) {
			std::cerr << ">>> Error: Domain1D::get_by_idx <<< Idx: " << i << " is out of domain " << _name << " of size " << _no_pts << std::endl;
		};
		return _min + i * _delta;
	};

	/********************
	Check if point is in domain
	********************/

	bool Domain1D::check_in_domain(double x) const
	{
		if (x < _min || x > _max) { 
			return false; 
		} else {
			return true;
		};
	};

	/********************
	Resize
	********************/

	void Domain1D::resize_no_pts_at_fixed_endpoints(int no_pts) {
		_no_pts = no_pts;
		// New delta
		_delta = (_max-_min)/(_no_pts-1);
	};

	void Domain1D::resize_no_pts_at_fixed_spacing(int no_pts) {
		_no_pts = no_pts;
		// New max
		_max = _min + (_no_pts-1)*_delta;
	};

	/********************
	Surrounding idxs
	********************/

	int Domain1D::get_surrounding_idxs(double x) const
	{
		int i = (x - _min) / _delta;
		if (i==_no_pts-1) {i--;};
		return i;
	};

	double Domain1D::get_frac_between(double x) const
	{
		return get_frac_between(x,get_surrounding_idxs(x));
	};

	double Domain1D::get_frac_between(double x, int i) const
	{
		return (x - get_by_idx(i)) / _delta;
	};
 
	/********************
	Print grid range
	********************/

	void Domain1D::print_domain_range() const
	{
		std::cout << "Domain: " << _name << " min: " << _min << " max: " << _max << std::endl;
	};

	/********************
	Write grid into an ofstream
	********************/

	void Domain1D::write_domain(std::string fname) const {
		std::ofstream f;
		f.open(fname);
		for (int i=0; i<_no_pts; i++) {
			f << _min + i * _delta << "\n";
		};
		f.close();
	};
};