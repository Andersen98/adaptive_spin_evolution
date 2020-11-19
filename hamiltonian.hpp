#include <array>
#include <vector>
#include <algorithm>
#include <complex>
#include <utility>
#include <cmath>
#include <iterator>


#include "input_tools.hpp"
#include "configuration.hpp"
#include <cassert>


//numeric
using std::complex;
using std::abs;

//containers
using std::pair;
using std::vector;
using std::array;

//algorithms
using std::sort;
using std::inplace_merge;
using std::set_difference;
using std::set_intersection;
using std::for_each;
using std::copy;


template<int num_modes, int num_bits>
class hamiltonian{

  typedef  State_Ket<bool,std::complex<double>,num_modes,num_bits> state_ket;
  typedef vector<state_ket> state_vector;


  array<array<double,num_modes>, (1<<num_bits)> g;
  array<double,num_modes> m;

  param_vals params;
  vector<int> mode_cap_exceeded;
  const int num_levels;


  state_vector psi_delta; //non-sorted (new terms added each step)
  state_vector psi_lbl; //label sorted
  state_vector psi_amp; //amp sorted
  
public:
  
  hamiltonian(state_vector initial_state ,
	      const param_vals &params_):g{{0}},m{0},params(params_),num_levels(1<<num_bits),psi_lbl(initial_state),psi_amp(initial_state){

    
    vector<pair<int,double>> idx_g_pairs(num_modes);
    for(int i=0;i< num_modes;i++){
      idx_g_pairs[i] = std::make_pair(i,params.mode_couplings[i]);
    }
    
    sort(idx_g_pairs.begin(),idx_g_pairs.end(),[](auto &it1,auto &it2){return it1.second > it2.second;});

    for(int j = 0; j < num_modes; j++){
	int i = idx_g_pairs[j].first;
	double gi = idx_g_pairs[j].second;
	//DEBUG
	assert(abs( gi - params.mode_couplings[i]) < .0001);
	//DEnd
	//descending order
	m[j] = params.mode_energies[i];
	
	for(int n = 0; n < num_levels; n++){
	  g[n][j] = gi*std::pow(n,.5);
	}
	//DEBUG
	assert(abs( g[1][j] - gi) < .0001);
	//DEnd
    }

    
  }
  
  
  static double abs_sqrd(const state_ket p){
    return(std::pow(std::real(p.amp),2)+std::pow(std::imag(p.amp),2));
  }
  static void normalize_state(state_vector &p){

    double N = 0;
    for(const state_ket &k: p){
      N += abs_sqrd(k);
    }
    N = 1/N;
    for_each(p.begin(),p.end(),[N](state_ket &k){k.amp*=N;});
    
  }

  state_vector get_psi()const{
    return psi_lbl;
  }
  int get_psi_size()const{
    return psi_lbl.size();
  }
    
  //returns current spin up, spin down 
  pair<double,double> get_spin_pop()const{
    pair<double,double> result = std::make_pair(0,0);

    for(const auto&k:psi_lbl){
      if(k.spin){
	result.first += abs_sqrd(k);
      }else{
	result.second += abs_sqrd(k);
      }	
    }
    
    return result;

  }

  //returns true to halt evaluation
  bool h_dipole(int idx,double dt){
    bool stop = false;

    double premagnitude = dt*abs_sqrd(psi_amp[idx]);
    if(g[num_levels-1][0]*premagnitude < params.energy_cutoff){
      return true;
    }

    //diagonal term
    // Sum omega_j n_j + params.emitter_energy;
    double sum_field=0;
    double sum_emitter = psi_amp[idx].spin ? params.up_energy:params.down_energy;
    int j = 0;
    complex<double> prefactor = -1i *dt * psi_amp[idx].amp;
    while((j < num_modes) && (params.energy_cutoff < premagnitude*g[num_levels-1][j]) ){

      int level = psi_amp[idx].get_mode(j);
      sum_field += level*m[j];
      complex<double> magnitude=0;
      switch(level){
      case(0):
	//ground state, only raise.
	if( (params.energy_cutoff < premagnitude*g[level+1][j] ) ){
	  magnitude = prefactor*g[level+1][j];
	  psi_delta.push_back(psi_amp[idx]);
	  psi_delta.back().spin = !psi_delta.back().spin;
	  psi_delta.back().amp = magnitude;
	  psi_delta.back().increment_mode(j);
	}
	break;
      case((1<<num_bits)-1):
	//max state, only lower
	if(params.energy_cutoff < premagnitude*g[level][j]){
	  magnitude = prefactor*g[level][j];
	  psi_delta.push_back(psi_amp[idx]);
	  psi_delta.back().spin = !psi_delta.back().spin;
	  psi_delta.back().amp = magnitude;
	  psi_delta.back().decrement_mode(j);
	}
	break;
      default:

	//raise
	if(params.energy_cutoff < premagnitude*g[level+1][j]){
	  magnitude = prefactor*g[level+1][j];
	  psi_delta.push_back(psi_amp[idx]);
	  psi_delta.back().spin = !psi_delta.back().spin;
	  psi_delta.back().amp = magnitude;
	  psi_delta.back().increment_mode(j);
	}
	//lower
	if(params.energy_cutoff < premagnitude*g[level][j]){
	  magnitude = prefactor*g[level][j];
	  psi_delta.push_back(psi_amp[idx]);
	  psi_delta.back().spin = !psi_delta.back().spin;
	  psi_delta.back().amp = magnitude;
	  psi_delta.back().decrement_mode(j);
	}
	break;


      }
      j++;
    }
    //compute diagonal term
    psi_amp[idx].amp += prefactor*(sum_field + sum_emitter);

    return(stop);

  }

  
  void merge_and_sort(){
    //psi_amp is the diagonal evolution
    //psi_delta is off diagonal evolution, but in general contains duplicates
    sort(psi_amp.begin(),psi_amp.end());
    sort(psi_delta.begin(),psi_delta.end());
    for( int i = 0; i < psi_amp.size(); i++){
      psi_lbl[i].amp += psi_delta[i].amp;
      //DEBUG Begin
      assert(psi_amp[i].spin==psi_lbl[i].spin);
      assert(psi_amp.size()==psi_lbl.size());
      //Dend
    }

     
    typename state_vector::iterator set_int,diff_end;
    psi_lbl.reserve(psi_lbl.size()+psi_delta.size());
    diff_end = psi_lbl.end();
    /* ========set_difference=============
      "Copies the elements from the sorted range [first1, last1) which are 
      not found in the sorted range [first2, last2) to the range beginning at d_first."
      cppreference
      returns "iterator past the end of the constructed range."
     
      Method: use set difference and std::back_inserter to append elements from delta that are not 
      in lbl
    */
    set_difference(psi_delta.begin(),psi_delta.end(),psi_lbl.begin(),psi_lbl.end(),std::back_inserter(psi_lbl));
    
    
    /*=======set intersection
      "Constructs a sorted range beginning at d_first consisting of elements that are found in 
      both sorted ranges [first1, last1) and [first2, last2). If some element is found m times 
      in [first1, last1) and n times in [first2, last2), the first std::min(m, n) elements
      will be copied from the first range to the destination range. The order of equivalent 
      elements is preserved. The resulting range cannot overlap with either of the input ranges. 

      returns Iterator past the end of the constructed range. 

      Method: the intersection is copied to psi_amp. then we add the two up to the consrtucted range
     */
  
    set_int= set_intersection(psi_delta.begin(),psi_delta.end(),psi_lbl.begin(),diff_end,
			      psi_amp.begin());
    int end_size = std::distance(psi_amp.begin(),set_int);
    for(int i = 0; i < end_size; i++){
      psi_lbl[i].amp += psi_amp[i].amp;
      
    }


    
    /*=========in place merge=================
      "Merges two consecutive sorted ranges [first, middle) and [middle, last)
      into one sorted range [first, last).A sequence [first, last) is said to be sorted with
      respect to a comparator comp if for any iterator it pointing to the sequence and any 
      non-negative integer n such that it + n is a valid iterator pointing to an element of 
      the sequence, comp(*(it + n), *it) evaluates to false. "
      
      returns(none)
     */
    inplace_merge(psi_lbl.begin(),diff_end,psi_lbl.end());

    psi_amp.resize(psi_lbl.size());
    normalize_state(psi_lbl);
    copy(psi_lbl.begin(),psi_lbl.end(),psi_amp.begin());
    sort(psi_amp.begin(),psi_amp.end(), [](auto &it1,auto &it2){return(abs_sqrd(it1)>abs_sqrd(it2));});   

  }

  void run(double dt){
    mode_cap_exceeded.clear();
    psi_delta.clear();
    psi_delta.reserve(num_modes*psi_amp.size());

    bool stop = false;
    for(int i = 0; i < psi_amp.size(); i++){
      stop = h_dipole(i,dt);
      if(stop){
	break;
	
      }
    }

    merge_and_sort();
    
  }
    
    
};
  

  
  

  
  
