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
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "input_tools.hpp"
#include "configuration.hpp"



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


//ublas

using boost::numeric::ublas::symmetric_adaptor;
class hamiltonian{
  
  typedef State_Ket<NUM_MODES,NUM_BITS> state_ket;
  typedef vector<state_ket> state_vector;
  typedef state_vector::iterator state_vector_iterator;
  //read write elements via lower triangle
  typedef boost::numeric::ublas::matrix<complex<double>> base_matrix_type;
  typedef symmetric_adaptor<base_matrix_type,boost::numeric::ublas::lower> matrix_type;
  typedef boost::numeric::ublas::vector<complex<double>> blas_vec;
  typedef boost::numeric::ublas::banded_matrix<int> banded_matrix_type;
  struct dir_edge_mode{
    uint out_idx;
    int connection_mode;
    bool raised;
  };
  
  struct ket_pair{
    state_ket raised;
    state_ket lowered;
    ket_pair();
    ket_pair(const state_ket &k);
  };
  
  array<array<double,NUM_MODES>, (1<<NUM_BITS)> g;
  array<double,NUM_MODES> m;
  param_vals params;
  array<int,NUM_MODES> new2old;
  array<int,NUM_MODES> old2new;
  
  
  const int num_levels;
      
  state_vector psi_delta; //non-sorted (new terms added each step)
  state_vector psi_lbl; //label sorted
  state_vector psi_amp; //amp sorted

  vector<vector<dir_edge_mode>> state_connections;
  uint matrix_idx_size;
  base_matrix_type base_matrix;
  matrix_type connection_matrix;
  banded_matrix_type spin_up_matrix;
  banded_matrix_type spin_down_matrix;
  blas_vec u;
  blas_vec v;
  //grow_configuration_space.cpp
  bool grow_configuration_space(int idx);

  //evolve_space.cpp
  void evolve_space(double dt);

  //helper used by append_connections
  //append_connections.cpp
  int binary_search_state(const state_ket &k,const state_vector &psi_search)const;
  void append_connections(state_vector &psi_mixed);
  void update_connection_matrix();
  //merge_states.cpp
  void merge_states();

  //setup.cpp
  void setup_connections();
  
  //core.cpp
  ket_pair get_connected_states(const state_ket &k,const int mode);
  static void normalize_state(state_vector &p);
  void update_matrix_diag(const state_ket &k);
  void update_matrix_diag_only_spin(uint update_idx,bool spin);
  void update_matrix_diag_only_mode(uint update_idx,int mode, int level);
  void grow_matrix(uint new_size);
public:
  array<int,NUM_MODES> mode_cap_exceeded;

  int evolve_state;
  constexpr static int grow_state = -1;
  constexpr static int blas_state0 = 0;
  constexpr static int blas_state1 = 1;

  //setup.cpp
  hamiltonian(const param_vals &params_);
  hamiltonian(const std::string &json_arg);
  

  //runtime.cpp
  void run_grow_evolve(double dt);
  void run_grow();
  void run_evolve(double dt);
  void set_epsilon(double e);
  void set_zero_except_init();
  void run_step(complex<double> factor);
  void switch_evolve();
  void blas_evolve( complex<double> factor);
  std::pair<double,double> get_blas_spin_pop()const;
  //odeint_evolve
  void odeint_evolve(double t1,double dt);

  //par_runtime.cpp
  void par_test_one();
  void par_test_two();
  
  //output.cpp
  pair<double,double> get_emitter_cavity_prob()const;
  state_vector get_state()const;
  matrix_type &get_matrix();
  int get_psi_size()const;
  pair<double,double> get_spin_pop()const;
  array<tuple<int,double,double>,NUM_MODES> get_modeLbl_quanta_pop()const;
  std::array<int,NUM_MODES> get_new2old()const;
  std::array<int,NUM_MODES> get_old2new()const;
  


};


#endif


