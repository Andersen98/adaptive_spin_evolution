#ifndef ADAPTIVE_SPIN_HAMILTONIAN
#define ADAPTIVE_SPIN_HAMILTONIAN 1


#include <cassert>
#include <array>
#include <vector>
#include <algorithm>
#include <complex>
#include <utility>
#include <cmath>
#include <iterator>
#include <tuple>


#include "../io_tools/input_tools.hpp"
#include "../configuration.hpp"



//numeric
using std::complex;
using std::abs;

//containers
using std::pair;
using std::vector;
using std::array;
using std::tuple;


//algorithms
using std::sort;
using std::inplace_merge;
using std::set_difference;
using std::set_intersection;
using std::for_each;
using std::copy;

class hamiltonian{

  typedef State_Ket<NUM_MODES,NUM_BITS> state_ket;
  typedef vector<state_ket> state_vector;
  typedef state_vector::iterator state_vector_iterator;

  array<array<double,NUM_MODES>, (1<<NUM_BITS)> g;
  array<double,NUM_MODES> m;
  array<int,NUM_MODES> new2old;
  array<int,NUM_MODES> old2new;
  param_vals params;
  vector<int> mode_cap_exceeded;
  const int num_levels;
  const complex<double> one_i;
      
  state_vector psi_delta; //non-sorted (new terms added each step)
  state_vector psi_lbl; //label sorted
  state_vector psi_amp; //amp sorted

  //keeps track of index
  int next_idx;
  struct dir_edge_mode{
    int out_idx;
    int connection_mode;
    bool raised;
  };
    
  vector<vector<dir_edge_mode>> state_connections;

  struct ket_pair{
    state_ket raised;
    state_ket lowered;
    ket_pair(){}
    ket_pair(const state_ket &k){
      raised = k;
      lowered = k;
    } 
    
  };
  


  //grow_configuration_space.cpp
  bool grow_configuration_space(int idx);

  //evolve_space.cpp
  void evolve_space(double dt);

  //helper used by append_connections
  //append_connections.cpp
  int binary_search_lbl(const state_ket &k);
  void append_connections(state_vector &delta);

  //merge_states.cpp
  void merge_states();

  //setup.cpp
  void setup_connections();
  
  //core.cpp
  ket_pair get_connected_states(const state_ket &k,const int mode);
  static double abs_sqrd(const state_ket &p);
  static void normalize_state(state_vector &p);
  static double abs(const state_ket &p);
public:
  //core.cpp
  void run_step(double dt);

  
  //setup.cpp
  hamiltonian(const param_vals &params_);
  
  //output.cpp
  state_vector &get_state_vector()const;
  int get_psi_size()const;
  pair<double,double> get_spin_pop()const;
  array<tuple<int,double,double>,NUM_MODES> get_modeLbl_quanta_pop()const;
  std::array<int,NUM_MODES> get_new2old()const;
  std::array<int,NUM_MODES> get_old2new()const;
  


};


#endif


