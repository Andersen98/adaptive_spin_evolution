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
#include <set>

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
	  psi_delta.back().amp = 0;
	  psi_delta.back().increment_mode(j);
	}
	break;
      case((1<<num_bits)-1):
	//max state, only lower
	if(params.energy_cutoff < premagnitude*g[level][j]){
	  psi_delta.push_back(psi_amp[idx]);
	  psi_delta.back().spin = !psi_delta.back().spin;
	  psi_delta.back().amp = 0;
	  psi_delta.back().decrement_mode(j);
	}
	break;
      default:

	//raise
	if(params.energy_cutoff < premagnitude*g[level+1][j]){
	  psi_delta.push_back(psi_amp[idx]);
	  psi_delta.back().spin = !psi_delta.back().spin;
	  psi_delta.back().amp = 0;
	  psi_delta.back().increment_mode(j);
	}
	//lower
	if(params.energy_cutoff < premagnitude*g[level][j]){
	  psi_delta.push_back(psi_amp[idx]);
	  psi_delta.back().spin = !psi_delta.back().spin;
	  psi_delta.back().amp = 0;
	  psi_delta.back().decrement_mode(j);
	}
	break;


      }
      j++;
    }
   
    return(stop);

  }

  void evolve_current_space(double dt){
    vector<complex<double>> delta(psi_amp.size(),complex<double>(0,0));
    std::pair<state_vector_iterator,state_vector_iterator> psi_el;
    for(int i = 0; i < psi_amp.size(); i++){
      complex<double> prefactor = complex<double>(0,-1)*dt*psi_amp[i].amp;
    
      
      for(int j = 0; j < num_modes; j++){
	
	int level = psi_amp[i].get_mode(j);
	
	//~~~~~Main std algorithm used ~~~~~~~~~~
	//For each mode in each ket, I attempt to evolve it in time within the current configuration
	//space. I think this could be optimized furthur, but right now I just want to see if it is sensible.
	
	//====lower bound (binary search)========
	/* ~~~~~~USAGE~~~~~~~~~~~~
	  ForwardIt lower_bound( ForwardIt first, ForwardIt last, const T& value );
	  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  ~~~~~~~~~~~DESCRIPTION~~~~~~~~~~~~~~~~~~
	  Returns an iterator pointing to the first element in the range [first, last) 
	  that is not less than (i.e. greater or equal to) value, or last if no such 
	  element is found.
	  
	  RETURNS: Iterator pointing to the first element that is not less than value, 
	  or last if no such element is found. 
	*/

	//~~~~~~Copy current ket to new ket
	state_ket new_ket = psi_amp[i];
	
	//~~~~~~Diagonal term E_atom + E_field~~~~~~~~~
	delta[i] += prefactor*(new_ket.spin ? params.up_energy:params.down_energy);
	delta[i] += prefactor*m[j]*double(level);
     
	//find if psi_amp has new_ket (first raised then lowered)
	switch(level){
	case(0):
	  //~~~~~~raise ket~~~~~~
	  new_ket.increment_mode(j);
	  new_ket.spin = !new_ket.spin;
	  psi_el = std::equal_range(psi_amp.begin(),psi_amp.end(), new_ket,
				    [](auto &el, auto &val){return el < val; });

	  if(psi_el.first != psi_el.second){
	    assert(std::distance(psi_el.first,psi_el.second) ==1);
	    int idx = std::distance(psi_amp.begin(),psi_el.first);
	    delta[idx] += prefactor*g[level+1][j];
	  }
	  break;
	case( (1<<num_bits) -1):
	  //~~~~~~~lower ket~~~~~~~~
	  new_ket.decrement_mode(j);
	  new_ket.spin = !new_ket.spin;
	  psi_el = std::equal_range(psi_amp.begin(),psi_amp.end(),new_ket,
				    [](auto &el, auto &val){return el < val; });
	  if(psi_el.first != psi_el.second){
	    int idx = std::distance(psi_amp.begin(),psi_el.first);
	    delta[idx] += prefactor*g[level][j];
	  }
	  break;  
	  
	default:
	  //~~~~~~~raise ket~~~~~~~~~~
	  new_ket.increment_mode(j);
	  new_ket.spin = !new_ket.spin;
	  psi_el = std::equal_range(psi_amp.begin(),psi_amp.end(),new_ket,
				    [](auto &el, auto &val){return el < val; });
	  if(psi_el.first != psi_el.second){
	    int idx = std::distance(psi_amp.begin(),psi_el.first);
	    delta[idx] += prefactor*g[level+1][j];
	  }
	  //~~~~~~~~lower ket~~~~~~~~~~~
	  new_ket.decrement_mode(j);
	  new_ket.decrement_mode(j);
     
	  psi_el = std::equal_range(psi_amp.begin(),psi_amp.end(),new_ket,
				    [](auto &el, auto &val){return el < val; });
	  if(psi_el.first != psi_el.second){
	    int idx = std::distance(psi_amp.begin(),psi_el.first);
	    delta[idx] += prefactor*g[level][j];
	  }
	  break;
	}//end switch
      
      }//end mode loop
    }//ket loop

    //add delta to amp
    for(int i =0; i < psi_amp.size(); i++){
      psi_amp[i].amp += delta[i];
    }
    

  }
    
  void union_and_evolve(double dt){
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
    psi_amp.resize(psi_delta.size()+psi_lbl.size());
    typename std::set<state_ket> s(psi_delta.begin(),psi_delta.end());
    psi_delta.assign(s.begin(),s.end());
    sort(psi_delta.begin(),psi_delta.end(),[](auto &it1, auto &it2){return it1 < it2; });
    
    psi_el = std::set_union(psi_lbl.begin(),psi_lbl.end(),psi_delta.begin(),psi_delta.end(),
			    psi_amp.begin(),[](auto &it1, auto &it2){return it1 < it2; });
    
    int amp_size = std::distance(psi_amp.begin(),psi_el);
    psi_amp.resize(amp_size);

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

    bool stop = false;
    for(int i = 0; i < psi_amp.size(); i++){
      stop = grow_configuration_space(i);
      if(stop){
	break;
	
      }
    }

    union_and_evolve(dt);
    
  }  
    
};
  

#endif  
