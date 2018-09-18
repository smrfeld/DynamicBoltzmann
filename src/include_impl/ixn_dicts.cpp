#include "../../include/dynamicboltz_bits/ixn_dicts.hpp"

// Other headers
#include "../../include/dynamicboltz_bits/species.hpp"
#include "../../include/dynamicboltz_bits/ixn_param.hpp"

#include <iostream>

/************************************
* Namespace for dblz
************************************/

namespace dblz {

	/****************************************
	BiasDict
	****************************************/

	// Constructor
	BiasDict::BiasDict() {};
	BiasDict::BiasDict(const BiasDict& other) {
		_copy(other);
	};
	BiasDict::BiasDict(BiasDict&& other) {
		_move(other);
	};
    BiasDict& BiasDict::operator=(const BiasDict& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    BiasDict& BiasDict::operator=(BiasDict&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };
	BiasDict::~BiasDict()
	{
		_clean_up();
	};
	void BiasDict::_clean_up() {};
	void BiasDict::_copy(const BiasDict& other) {
		_dict = other._dict;
	};
	void BiasDict::_move(BiasDict& other) {
		_dict = other._dict;
		other._dict.clear();
	};

	// Add to the dict
	void BiasDict::add_ixn(Sptr sp, Iptr ixn) {
		_dict[sp].push_back(ixn);
	};

	// Get from the dict
	double BiasDict::get_ixn_at_timepoint(Sptr sp, int timepoint) const {
		auto it1 = _dict.find(sp);
		if (it1 == _dict.end()) {
			return 0.0;
		};

		double r=0.0;
		for (auto ptr: it1->second) {
			r += ptr->get_val_at_timepoint(timepoint);
		};
		return r;
	};






























	/****************************************
	O2IxnDict
	****************************************/

	class O2IxnDict::Impl {

	private:

		// Site to site to ixns
		std::map<Sptr,std::map<Sptr,std::vector<Iptr>>> _dict;

		// Add ixn
		void _add_ixn(Sptr sp_of_site_1, Sptr sp_of_site_2, Iptr ixn);

		// Constructor helpers
		void _clean_up();
		void _copy(const Impl& other);
		void _move(Impl& other);

	public:

		// Constructor		
		Impl();
		Impl(const Impl& other);
		Impl(Impl&& other);
		Impl& operator=(const Impl& other);
		Impl& operator=(Impl&& other);
		~Impl();

		// Add to the dict
		void add_ixn(Sptr sp_of_site_1, Sptr sp_of_site_2, Iptr ixn, bool this_order=false);

		// Get from the dict
		double get_ixn_at_timepoint(Sptr sp_of_site_1, Sptr sp_of_site_2, int timepoint) const;
	};

	// Constructor
	O2IxnDict::Impl::Impl() {};
	O2IxnDict::Impl::Impl(const Impl& other) {
		_copy(other);
	};
	O2IxnDict::Impl::Impl(Impl&& other) {
		_move(other);
	};
    O2IxnDict::Impl& O2IxnDict::Impl::operator=(const Impl& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    O2IxnDict::Impl& O2IxnDict::Impl::operator=(Impl&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };
	O2IxnDict::Impl::~Impl()
	{
		_clean_up();
	};
	void O2IxnDict::Impl::_clean_up() {};
	void O2IxnDict::Impl::_copy(const Impl& other) {
		_dict = other._dict;
	};
	void O2IxnDict::Impl::_move(Impl& other) {
		_dict = other._dict;
		other._dict.clear();
	};

	// Add to the dict
	void O2IxnDict::Impl::_add_ixn(Sptr sp_of_site_1, Sptr sp_of_site_2, Iptr ixn) {
		auto it = std::find(_dict[sp_of_site_1][sp_of_site_2].begin(),_dict[sp_of_site_1][sp_of_site_2].end(),ixn);
		if (it == _dict[sp_of_site_1][sp_of_site_2].end()) {
			_dict[sp_of_site_1][sp_of_site_2].push_back(ixn);
		};
	};
	void O2IxnDict::Impl::add_ixn(Sptr sp1, Sptr sp2, Iptr ixn, bool this_order) {
		_add_ixn(sp1,sp2,ixn);
		if (!this_order) {
			_add_ixn(sp2,sp1,ixn);
		};
	};

	// Get from the dict
	double O2IxnDict::Impl::get_ixn_at_timepoint(Sptr sp_of_site_1, Sptr sp_of_site_2, int timepoint) const {
		auto it1 = _dict.find(sp_of_site_1);
		if (it1 == _dict.end()) {
			return 0.0;
		};
		auto it2 = it1->second.find(sp_of_site_2);
		if (it2 == it1->second.end()) {
			return 0.0;
		};

		double r=0.0;
		for (auto ptr: it2->second) {
			r += ptr->get_val_at_timepoint(timepoint);
		};
		return r;
	};

	// Constructor forward
	O2IxnDict::O2IxnDict() : _impl(new Impl()) {};
	O2IxnDict::O2IxnDict(const O2IxnDict& other) : _impl(new Impl(*other._impl)) {};
	O2IxnDict::O2IxnDict(O2IxnDict&& other) : _impl(std::move(other._impl)) {};
	O2IxnDict& O2IxnDict::operator=(const O2IxnDict& other) {
		_impl.reset(new Impl(*other._impl));
		return *this;
	};
	O2IxnDict& O2IxnDict::operator=(O2IxnDict&& other) {
		_impl = std::move(other._impl);
		return *this;
	};
	O2IxnDict::~O2IxnDict() = default;

	// Add to the dict forward
	void O2IxnDict::add_ixn(Sptr sp_of_site_1, Sptr sp_of_site_2, Iptr ixn, bool this_order) {
		_impl->add_ixn(sp_of_site_1,sp_of_site_2,ixn,this_order);
	};

	// Get from the dict forward
	double O2IxnDict::get_ixn_at_timepoint(Sptr sp_of_site_1, Sptr sp_of_site_2, int timepoint) const {
		return _impl->get_ixn_at_timepoint(sp_of_site_1,sp_of_site_2,timepoint);
	};
























	/****************************************
	O3IxnDict
	****************************************/

	class O3IxnDict::Impl {

	private:

		// Site 1 to site 2 to site 3 to ixns
		std::map<Sptr,std::map<Sptr,std::map<Sptr,std::vector<Iptr>>>> _dict;

		// Add to the dict internal
		void _add_ixn(Sptr sp_of_site_1, Sptr sp_of_site_2, Sptr sp_of_site_3, Iptr ixn);

		// Constructor helpers
		void _clean_up();
		void _copy(const Impl& other);
		void _move(Impl& other);

	public:

		// Constructor		
		Impl();
		Impl(const Impl& other);
		Impl(Impl&& other);
		Impl& operator=(const Impl& other);
		Impl& operator=(Impl&& other);
		~Impl();

		// Add to the dict
		void add_ixn(Sptr sp_of_site_1, Sptr sp_of_site_2, Sptr sp_of_site_3, Iptr ixn, bool this_order);

		// Get from the dict
		double get_ixn_at_timepoint(Sptr sp_of_site_1, Sptr sp_of_site_2, Sptr sp_of_site_3, int timepoint) const;
	};

	// Constructor
	O3IxnDict::Impl::Impl() {};
	O3IxnDict::Impl::Impl(const Impl& other) {
		_copy(other);
	};
	O3IxnDict::Impl::Impl(Impl&& other) {
		_move(other);
	};
    O3IxnDict::Impl& O3IxnDict::Impl::operator=(const Impl& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    O3IxnDict::Impl& O3IxnDict::Impl::operator=(Impl&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };
	O3IxnDict::Impl::~Impl()
	{
		_clean_up();
	};
	void O3IxnDict::Impl::_clean_up() {};
	void O3IxnDict::Impl::_copy(const Impl& other) {
		_dict = other._dict;
	};
	void O3IxnDict::Impl::_move(Impl& other) {
		_dict = other._dict;
		other._dict.clear();
	};

	// Add to the dict
	void O3IxnDict::Impl::_add_ixn(Sptr sp_of_site_1, Sptr sp_of_site_2, Sptr sp_of_site_3, Iptr ixn) {
		auto it = std::find(_dict[sp_of_site_1][sp_of_site_2][sp_of_site_3].begin(),_dict[sp_of_site_1][sp_of_site_2][sp_of_site_3].end(),ixn);
		if (it == _dict[sp_of_site_1][sp_of_site_2][sp_of_site_3].end()) {
			_dict[sp_of_site_1][sp_of_site_2][sp_of_site_3].push_back(ixn);
		};
	};
	void O3IxnDict::Impl::add_ixn(Sptr sp1, Sptr sp2, Sptr sp3, Iptr ixn, bool this_order) {
		_add_ixn(sp1,sp2,sp3,ixn);
		if (!this_order) {
			_add_ixn(sp3,sp2,sp1,ixn);
		};
	};

	// Get from the dict
	double O3IxnDict::Impl::get_ixn_at_timepoint(Sptr sp_of_site_1, Sptr sp_of_site_2, Sptr sp_of_site_3, int timepoint) const {
		auto it1 = _dict.find(sp_of_site_1);
		if (it1 == _dict.end()) {
			return 0.0;
		};
		auto it2 = it1->second.find(sp_of_site_2);
		if (it2 == it1->second.end()) {
			return 0.0;
		};
		auto it3 = it2->second.find(sp_of_site_3);
		if (it3 == it2->second.end()) {
			return 0.0;
		};

		double r=0.0;
		for (auto ptr: it3->second) {
			r += ptr->get_val_at_timepoint(timepoint);
		};
		return r;
	};

	// Constructor forward
	O3IxnDict::O3IxnDict() : _impl(new Impl()) {};
	O3IxnDict::O3IxnDict(const O3IxnDict& other) : _impl(new Impl(*other._impl)) {};
	O3IxnDict::O3IxnDict(O3IxnDict&& other) : _impl(std::move(other._impl)) {};
	O3IxnDict& O3IxnDict::operator=(const O3IxnDict& other) {
		_impl.reset(new Impl(*other._impl));
		return *this;
	};
	O3IxnDict& O3IxnDict::operator=(O3IxnDict&& other) {
		_impl = std::move(other._impl);
		return *this;
	};
	O3IxnDict::~O3IxnDict() = default;

	// Add to the dict forward
	void O3IxnDict::add_ixn(Sptr sp_of_site_1, Sptr sp_of_site_2, Sptr sp_of_site_3, Iptr ixn, bool this_order) {
		_impl->add_ixn(sp_of_site_1,sp_of_site_2,sp_of_site_3,ixn,this_order);
	};

	// Get from the dict forward
	double O3IxnDict::get_ixn_at_timepoint(Sptr sp_of_site_1, Sptr sp_of_site_2, Sptr sp_of_site_3, int timepoint) const {
		return _impl->get_ixn_at_timepoint(sp_of_site_1,sp_of_site_2,sp_of_site_3,timepoint);
	};
};