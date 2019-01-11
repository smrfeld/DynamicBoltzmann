
.. _program_listing_file_lattice.hpp:

Program Listing for File lattice.hpp
====================================

- Return to documentation for :ref:`file_lattice.hpp`

.. code-block:: cpp

   // string
   #ifndef STRING_h
   #define STRING_h
   #include <string>
   #endif
   
   // vector
   #ifndef VECTOR_h
   #define VECTOR_h
   #include <vector>
   #endif
   
   // list
   #ifndef LIST_h
   #define LIST_h
   #include <list>
   #endif
   
   // map
   #ifndef MAP_h
   #define MAP_h
   #include <map>
   #endif
   
   /************************************
   * Namespace for DynamicBoltzmann
   ************************************/
   
   namespace DynamicBoltzmann {
   
       /****************************************
       General functions
       ****************************************/
   
       struct Species;
       struct Site;
       typedef std::list<Site> lattice;
       typedef std::list<Site>::iterator latt_it;
   
       /****************************************
       Structure to hold a lattice site
       ****************************************/
   
       struct Site {
           int x;
           int y;
           int z;  
           Species *sp;
           std::vector<latt_it> nbrs;
   
           // Constructor
           Site();
           Site(int xIn, int yIn, int zIn);
           Site(int xIn, int yIn, int zIn, Species *spIn);
       };
       // Comparator
       bool operator <(const Site& a, const Site& b);
       bool operator==(const Site& a, const Site& b);
       std::ostream& operator<<(std::ostream& os, const Site& s);
   
       /****************************************
       Lattice
       ****************************************/
   
       class Lattice
       {
       private:
   
           // Internal maps
           lattice _latt;
   
           // Size
           int _box_length;
   
           // Pointers to species present
           std::map<std::string,Species*> _sp_map;
   
           /********************
           Lookup a site iterator from x,y,z
           ********************/
   
           latt_it _look_up(int x, int y, int z);
   
       public:
   
           /********************
           Constructor/Destructor
           ********************/
   
           Lattice(int box_length);
           Lattice();
           ~Lattice();
   
           /********************
           Add a species
           ********************/
   
           void add_species(Species *sp);
   
           /********************
           Clear, size
           ********************/
   
           void clear();
           int size();
   
           /********************
           Make a mol
           ********************/
   
           bool make_mol(latt_it s, Species *sp);
   
           /********************
           Erase a mol
           ********************/
   
           bool erase_mol(latt_it s);
   
           /********************
           Write lattice to a file
           ********************/
   
           void write_to_file(std::string fname);
   
           /********************
           Read lattice from a file
           ********************/
   
           void read_from_file(std::string fname);
   
           /********************
           Anneal
           ********************/
   
           void anneal(int n_steps);
       };
   
   };
