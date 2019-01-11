
.. _program_listing_file_species.cpp:

Program Listing for File species.cpp
====================================

- Return to documentation for :ref:`file_species.cpp`

.. code-block:: cpp

   #include <iostream>
   #include "ixn_param.hpp" // also includes species header
   
   /************************************
   * Namespace for Gillespie3D
   ************************************/
   
   namespace DynamicBoltzmann {
   
       /****************************************
       Species
       ****************************************/
       
       Species::Species(std::string name) {
           _name = name;
           _count = 0;
   
           // Ptrs
           _t_opt_ptr = nullptr;
           _h_ptr = nullptr;
       };
   
       /********************
       Set pointer to the opt time variable
       ********************/
   
       void Species::set_opt_time_ptr(int *t_opt_ptr)
       {
           _t_opt_ptr = t_opt_ptr;
       };
   
       /********************
       Set h, j ptr
       ********************/
   
       void Species::set_h_ptr(IxnParam *h_ptr) {
           _h_ptr = h_ptr;
       };
       void Species::add_j_ptr(Species* sp, IxnParam *j_ptr) {
           _j_ptr[sp] = j_ptr;
           // Also add entry in nn count
           _nn_count[sp] = 0;
       };
   
       /********************
       Validate setup
       ********************/
   
       void Species::validate_setup() const {
           std::cout << "--- Validate species: " << _name << " ---" << std::endl;
           std::cout << "   NNs: " << std::flush;
           for (auto p: _nn_count) {
               std::cout << p.first->name() << " " << std::flush;
           };
           std::cout << std::endl;
           if (!_t_opt_ptr) {
               std::cerr << "ERROR: no time ptr set" << std::endl;
               exit(EXIT_FAILURE);
           } else {
               std::cout << "   Time ptr is set" << std::endl;
           };
           if (!_h_ptr) {
               std::cerr << "ERROR: no h ptr set" << std::endl;
               exit(EXIT_FAILURE);
           } else {
               std::cout << "   h ptr is set" << std::endl;
           };
       };
   
       /********************
       Setters/getters
       ********************/
   
       double Species::h() const {
           return _h_ptr->get_at_time(*_t_opt_ptr);
       };
       double Species::j(Species *other) const {
           return _j_ptr.at(other)->get_at_time(*_t_opt_ptr);
       };
       int Species::count() const {
           return _count;
       };
       int Species::nn_count(Species *other) const {
           return _nn_count.at(other);
       };
   
       std::string Species::name() const { return _name; };
   
       /********************
       Increment counts
       ********************/
   
       void Species::count_plus() { _count++; };
       void Species::count_minus() { _count--; };
       void Species::nn_count_plus(Species* other) { _nn_count[other]++; };
       void Species::nn_count_minus(Species* other) { _nn_count[other]--; };
   
       /********************
       Reset counts
       ********************/
   
       void Species::reset_counts() {
           _count = 0;
           for (auto it=_nn_count.begin(); it!=_nn_count.end(); it++) {
               it->second = 0;
           };
       };
   
       /********************
       Comparator
       ********************/
   
       bool operator <(const Species& a, const Species& b) {
           return a.name() < b.name();
       };
   
   };
