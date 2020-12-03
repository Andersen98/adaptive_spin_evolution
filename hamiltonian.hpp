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


template<int num_modes, int num_bits>
class hamiltonian{

  typedef State_Ket<bool,std::complex<double>,num_modes,num_bits> state_ket;
  typedef vector<state_ket> state_vector;
  using state_vector_iterator = typename state_vector::iterator;

  array<array<double,num_modes>, (1<<num_bits)> g;
  array<double,num_modes> m;

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
  
  
public:
  
  hamiltonian(const param_vals &params_):one_i(0,1),g{{0}},m{0},params(params_),num_levels(1<<num_bits),psi_lbl(0),psi_amp(0),next_idx(0),state_connections(){

       
    vector<pair<int,double>> idx_g_pairs(num_modes);
    vector<int> old2new (num_modes);
    for(int i=0;i< num_modes;i++){
      idx_g_pairs[i] = std::make_pair(i,params.mode_couplings[i]);
    }
    
    sort(idx_g_pairs.begin(),idx_g_pairs.end(),[](auto &it1,auto &it2){return it1.second > it2.second;});

    //map from old idx to new idx
    for(int i = 0; i < num_modes; i++){
      int j =  idx_g_pairs[i].first;
      old2new[j] = i;
    }

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

    //construct initial states
    for(auto &sk: params.initial_state){
      int mode = old2new[sk.mode];
      state_ket k;
      k.set_mode(mode,sk.n);
      k.spin = sk.spin;
      k.amp = sk.amp;
      psi_delta.push_back(k);
    }

    //start with a single root state.
    sort(psi_delta.begin(),psi_delta.end(),[](auto &it1,auto &it2){return it1<it2;});
    psi_lbl.push_back(psi_delta[0]);
    psi_lbl[0].idx = next_idx++;

    
    //do a run for zero seconds
    union_and_evolve(0);
    
    normalize_state(psi_lbl);
    sort(psi_lbl.begin(),psi_lbl.end(),[](auto &it1,auto &it2){return it1<it2;});
    psi_amp.resize(psi_lbl.size());
    copy(psi_lbl.begin(),psi_lbl.end(),psi_amp.begin());
    sort(psi_amp.begin(),psi_amp.end(),[](auto &it1,auto &it2){return abs_sqrd(it1) > abs_sqrd(it2);});
    
    
  }
  
  
  static double abs_sqrd(const state_ket p){
    return(std::pow(std::real(p.amp),2)+std::pow(std::imag(p.amp),2));
  }
  static void normalize_state(state_vector &p){

    double N = 0;
    for(const state_ket &k: p){
      N += abs_sqrd(k);
    }
    N = 1/std::sqrt(N);
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
  
  array<tuple<int,double,double>,num_modes> get_modeLbl_quanta_pop(){
    array<tuple<int,double,double>,num_modes> result;

    for_each(result.begin(),result.end(),[j=0](auto &el)mutable{el= std::make_tuple<int,double,double>(j++,0,0);});
    for_each(result.begin(),result.end(),[&](auto &el){
      for(auto &ket:psi_lbl){
	double probability =abs_sqrd(ket);
	int number_expectation =0;
	int modeLbl = std::get<0>(el);
	int nVal = ket.get_mode(modeLbl);
	std::get<1>(el) += nVal*probability;
	std::get<2>(el) += probability;
	  
      };});
    return result;
  }

         
    
  //returns true to halt evaluation
  bool grow_configuration_space(int idx){
    bool stop = false;

    double premagnitude = abs_sqrd(psi_amp[idx]);
    if(g[num_levels-1][0]*premagnitude < params.energy_cutoff){
      return true;
    }

    //diagonal term
    // Sum omega_j n_j + params.emitter_energy;
    int j = 0;
    complex<double> prefactor = std::complex<double>(0,-1) * psi_amp[idx].amp;
    while((j < num_modes) && (params.energy_cutoff < premagnitude*g[num_levels-1][j]) ){

      int level = psi_amp[idx].get_mode(j);
      switch(level){
      case(0):
	//ground state, only raise.
	if( (params.energy_cutoff < premagnitude*g[level+1][j] ) ){
	  psi_delta.push_back(psi_amp[idx]);
	  psi_delta.back().spin = !psi_delta.back().spin;
	  psi_delta.back().increment_mode(j);
	}
	break;
      case((1<<num_bits)-1):
	//max state, only lower
	if(params.energy_cutoff < premagnitude*g[level][j]){
	  psi_delta.push_back(psi_amp[idx]);
	  psi_delta.back().spin = !psi_delta.back().spin;
	  psi_delta.back().decrement_mode(j);
	}
	break;
      default:
	
	bool should_raise = params.energy_cutoff < premagnitude*g[level+1][j];
	bool should_lower = params.energy_cutoff < premagnitude*g[level][j];
	
	if(should_raise && should_lower ){
	  //raise
	  psi_delta.push_back(psi_amp[idx]);
	  psi_delta.back().spin = !psi_delta.back().spin;
	  psi_delta.back().increment_mode(j);
	  //lower
	  psi_delta.push_back(psi_amp[idx]);
	  psi_delta.back().spin = !psi_delta.back().spin;
	  psi_delta.back().decrement_mode(j);
	}
	break;


      }
      j++;
    }
   
    return(stop);

  }

  
  void evolve_current_space(double dt){

    //sort by label

    vector<int> con2amp(psi_amp.size());
    for(int i = 0; i < psi_amp.size(); i++){
      //con -> ampidx
      con2amp[psi_amp[i].idx] = i;
    }
    vector<complex<double>> delta(psi_amp.size());

    for(int i = 0; i < next_idx; i ++){
      complex<double> start_amp = psi_amp[con2amp[i]].amp;
      //diagonal term
      //--atom energy
      delta[i] += start_amp*complex<double>(0,-1)*dt*(psi_amp[con2amp[i]].spin? params.up_energy:params.down_energy);

      assert(i==psi_amp[con2amp[i]].idx);
      for(auto & edge:state_connections[i]){
	int out_idx = edge.out_idx;
	complex<double> out_amp = psi_amp[con2amp[out_idx]].amp;
	int connection_mode = edge.connection_mode;
	int start_level = psi_amp[con2amp[i]].get_mode(connection_mode);
	bool raised =  edge.raised;
	//start -> finish
	double g_factor = raised?(g[start_level+1][connection_mode]):(g[start_level][connection_mode]);
	delta[out_idx] += start_amp*complex<double>(0,-1)*dt*g_factor;
	//diagonal term
	delta[i] += start_amp*complex<double>(0,-1)*dt*double(start_level)*m[connection_mode];

	//finish->start
	delta[i] += out_amp*complex<double>(0,-1)*dt*g_factor;

      }
    }


    for(int i = 0; i < psi_amp.size(); i++){
      psi_amp[con2amp[i]].amp += delta[i];
      
    }

  }


  ket_pair get_connected_states(const state_ket &k,const int mode){

    int level = k.get_mode(mode);
    ket_pair kp(k);
    kp.raised.idx = state_ket::null_idx;
    kp.lowered.idx = state_ket::null_idx;
    
    switch(level){
     case(0):
       //raise
       kp.raised.idx = state_ket::empty_idx;
       kp.raised.increment_mode(mode);
       kp.raised.spin = !kp.raised.spin;
  
       break;
     case((1<<num_bits)-1):
       //lower
       kp.lowered.idx = state_ket::empty_idx;
       kp.lowered.decrement_mode(mode);
       kp.lowered.spin = !kp.lowered.spin;
       
       break;
     default:
       //raise
       kp.raised.idx = state_ket::empty_idx;
       kp.raised.increment_mode(mode);
       kp.raised.spin = !kp.raised.spin;

       //lower
       kp.lowered.idx = state_ket::empty_idx;
       kp.lowered.decrement_mode(mode);
       kp.lowered.spin = !kp.lowered.spin;
       break;


     }//end switch

    return kp;

  }

  int binary_search_lbl(const state_ket &k){
    std::pair<state_vector_iterator,state_vector_iterator> psi_el;
    int result = -1; //if it fails, the result is -1;
    psi_el = std::equal_range(psi_lbl.begin(),psi_lbl.end(),k,
			      [](auto &el,auto &val){return el < val;});

    if(psi_el.first != psi_el.second){
      result = std::distance(psi_lbl.begin(),psi_el.first);
    }

    return(result);
  }

  void setup_connections(state_vector &psi){

    for(auto &k:psi){
      k.idx = next_idx++;
    }

     
    for(auto &k: psi){
      vector<dir_edge_mode> edges;
              
      //apply hamiltonian
      for(int i = 0; i < num_modes; i++){
	
	ket_pair kp = get_connected_states(k,i);
	bool raised = kp.raised.idx==state_ket::empty_idx;
	bool lowered = kp.lowered.idx==state_ket::empty_idx;
	 
	if(raised){
	  //look for an instance
	  int vec_idx = binary_search_lbl(kp.raised);
	  if(vec_idx>0){
	    dir_edge_mode em;
	    em.out_idx = psi_lbl[vec_idx].idx;
	    em.connection_mode = i;
	    em.raised = true;
	    edges.push_back(em);
	  }
	  
	}
	if(lowered){
	  int vec_idx = binary_search_lbl(kp.lowered);
	  if(vec_idx){
	    dir_edge_mode em;
	    em.out_idx = psi_lbl[vec_idx].idx;
	    em.connection_mode = i;
	    em.raised = false;
	    edges.push_back(em); 
	  }
	}
	
      }//end mode loop
      
      state_connections.push_back(edges);
      
    }//end k loop
    
  }
  void append_connections(state_vector &delta){

   
    
     for(auto &k: delta){
       vector<dir_edge_mode> edges;
       k.idx = next_idx++;
       
       //apply hamiltonian
       for(int i = 0; i < num_modes; i++){
	 
	 ket_pair kp = get_connected_states(k,i);
	 bool raised = kp.raised.idx==state_ket::empty_idx;
	 bool lowered = kp.lowered.idx==state_ket::empty_idx;

	 if(raised){
	   //look for an instance
	   int vec_idx = binary_search_lbl(kp.raised);
	   if(vec_idx>0){
	     dir_edge_mode em;
	     em.out_idx = psi_lbl[vec_idx].idx;
	     em.connection_mode = i;
	     em.raised = true;
	     edges.push_back(em);
	   }

	 }
	 if(lowered){
	   int vec_idx = binary_search_lbl(kp.lowered);
	   if(vec_idx){
	     dir_edge_mode em;
	     em.out_idx = psi_lbl[vec_idx].idx;
	     em.connection_mode = i;
	     em.raised = false;
	     edges.push_back(em); 
	   }
	 }

       }//end mode loop

       state_connections.push_back(edges);
       
     }//end k loop
  }
     


 
  void union_and_evolve(double dt){
    if(psi_delta.size()>0){
      //==============set_union====================
      /*USAGE
	template< class InputIt1, class InputIt2, class OutputIt >
	OutputIt set_union( InputIt1 first1, InputIt1 last1,InputIt2 first2, InputIt2 last2,
	                       OutputIt d_first );
	**DESCRIPTION
	Constructs a sorted union beginning at d_first consisting of the set of elements present 
	in one or both sorted ranges [first1, last1) and [first2, last2). 
	**RETURNS
	Iterator past the end of the constructed range. 
      */
      typename state_vector::iterator psi_el;

    
      sort(psi_delta.begin(),psi_delta.end(),[](auto &it1, auto &it2){return it1 < it2; });
      //dedupe
      vector<state_ket> delta_clean(psi_delta.size());
      delta_clean[0] = psi_delta[0];
      int j = 1;
      for(int i = 1; i < psi_delta.size(); i++){
	if(delta_clean[j-1] < psi_delta[i]){
	  delta_clean[j] = psi_delta[i];
	}	
      }
      //RESIZE DELTA CLEAN BEFORE DOING ANYTHING ELSE!
      int delta_clean_size = j;
      delta_clean.resize(delta_clean_size);

      
      //append new connections before merge
      append_connections(delta_clean);

      //merge into psi_lbl
      psi_amp.resize(delta_clean.size()+psi_lbl.size());    
      psi_el = std::set_union(psi_lbl.begin(),psi_lbl.end(),delta_clean.begin(),delta_clean.end(),
			      psi_amp.begin(),[](auto &it1, auto &it2){return it1 < it2; });
      
    }//endif delta_size.()
    
    //evolves psi_amp
    evolve_current_space(dt);
    //move psi_amp to psi_lbl
    normalize_state(psi_amp);
    psi_lbl.resize(psi_amp.size());
    std::copy(psi_amp.begin(),psi_amp.end(),psi_lbl.begin());
    std::sort(psi_amp.begin(),psi_amp.end(),[](state_ket&it1,state_ket&it2){return abs_sqrd(it1)>abs_sqrd(it2);});
    
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

  void run_step(double dt){
    mode_cap_exceeded.clear();
    psi_delta.clear();
    psi_delta.reserve(num_modes*psi_amp.size());

    /*
    bool stop = false;
    for(int i = 0; i < psi_amp.size(); i++){
      stop = grow_configuration_space(i);
      if(stop){
	break;
	
      }
    }
    */
    union_and_evolve(dt);
    
  }  
    
};
  

#endif  
