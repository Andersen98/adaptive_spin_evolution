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
#include <memory>


#include "input_tools.hpp"
#include "configuration.hpp"
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/SparseCore>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

//memory
using std::shared_ptr;

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
  typedef vector<shared_ptr<state_ket>> state_vector;
  typedef state_vector::iterator state_vector_iterator;
  
  //read write elements via lower triangle

  
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
  typedef Eigen::Triplet<std::complex<double>> T;
  typedef Eigen::SparseMatrix<std::complex<double>> SpMat;
  typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic,1> ComplexVec;
  typedef Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> ComplexMatrix;
  SpMat H_matrix;
  ComplexVec psi_u;
  ComplexVec psi_uinit;
  ComplexMatrix H_eigen_vectors;
  ComplexVec H_eigen_vals;
  
  
  Eigen::Matrix<std::complex<double>,2, Eigen::Dynamic> SpinMatrix;
  std::vector<T>  state_connections;
  int next_idx;
  
  //grow_configuration_space.cpp
  bool grow_configuration_space(int idx);

 
  void add_connection(state_ket &k, state_vector &psi_mixed);
  //merge_states.cpp
  void merge_states();

  //setup.cpp
  void setup_connections();
  
  //core.cpp
  ket_pair get_connected_states(const state_ket &k,const int mode);
  static void normalize_state(state_vector &p);
  std::complex<double> matrix_diag(const state_ket &k);
  std::complex<double> matrix_diag_only_spin(bool spin);
  std::complex<double> matrix_diag_only_mode(int mode, int level);

public:
  array<int,NUM_MODES> mode_cap_exceeded;

  //setup.cpp
  hamiltonian(const param_vals &params_);
  hamiltonian(const std::string &json_arg);
  
   //helper used by append_connections
  //append_connections.cpp
  int binary_search_state(const state_ket &k,const state_vector &psi_search)const;
  //searched psi_lbl for a particular ket
  int get_ket_idx(const state_ket &k) const;
  //runtime.cpp
  void reset();
  void reset_with_state(state_vector init_state);
  void run_grow();
  void set_epsilon(double e);
  void set_zero_except_init();
  void run_step(complex<double> factor);
  void store_matrix();
  void store_vector();
  std::vector<std::tuple<int,int,std::complex<double>>> get_tuples()const;
  std::vector<std::pair<int, int>> get_spin_idxs()const;
  pair<double,double> evolve_state(double time);
  pair<double,double> evolve_step(complex<double> factor);
  //par_runtime.cpp
  void par_test_one();
  void par_test_two();

  //output.cpp
  pair<double,double> get_emitter_cavity_prob()const;
  state_vector get_state()const;
  int get_psi_size()const;
  pair<double,double> get_spin_pop()const;
  array<tuple<int,double,double>,NUM_MODES> get_modeLbl_quanta_pop()const;
  std::array<int,NUM_MODES> get_new2old()const;
  std::array<int,NUM_MODES> get_old2new()const;
  


};


#endif


