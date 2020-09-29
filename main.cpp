#include "configuration.hpp"
#include "input_tools.hpp"
#include <vector>
#include <iostream>
#include <bitset>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <cmath>
#include <complex>
#include <utility>
#include <string>
#include <cassert>
using namespace std;


int main(int argc, char * argv[]){
  
  //parameters to store system information
  //parameters also store simulation parameters like
  //energy tolerance and time step duration
  param_vals params;

  try{
    bool exit = get_params(params,argc,argv);
    if(exit)
      return 0;
  }catch(std::exception& e){
    cout << e.what() << endl;
    return 1;
  }
    params.energy_cutoff;

  //make transitions
  vector<int> s_idx(params.mode_energies.size());
  iota(s_idx.begin(),s_idx.end(),0);

  sort(s_idx.begin(),s_idx.end(),
       [&x=params.mode_couplings](auto i1,auto i2){ return x[i1] > x[i2];});

  //make full transition list 
  vector<vector<double>> g(params.max_occupation);
  transform(g.begin(),g.end(),g.begin(),
	    [&m=params.mode_couplings, j=1](vector<double> el) mutable {
	      vector<double> v(m);
	      transform(v.begin(),v.end(),v.begin(), [&j](double g){return sqrt(j)*g;});
	      j++;
	      return v;
	    });
    

  const int epsilon = params.energy_cutoff;
  const int num_spins = 2;
  const int num_modes = 100;
  const int num_bits = 4;
  const int max_level = (1<<num_bits) -1;
  const int giant_words = (num_bits*num_modes)/64 +1;
  typedef partial_config<bool,num_modes,num_bits,giant_words> full_config;
  typedef pair<full_config,complex<double>> state_pair;
  typedef vector<state_pair> state_vec;
  state_vec psi_lbl(1);
  psi_lbl[0].first.set_spin(true);
  psi_lbl[0].second = 1;
  state_vec psi_amp(psi_lbl);
  double t = params.t0;
  double dt = (params.tf-params.t0)/( double(params.N));
  assert(g[0].size() == num_modes);
  while(t < params.tf){
    
    //iterate through each component of state vector
    //psi_n+1 = (1 - i H dt)psi_n
    //        = psi_n - i H dt psi_n
    //delta_psi = (-i H dt) psi_n
    //psi_n = \sum spin (x) [mode_levels]
    
    state_vec delta_psi;
    
    for(auto &conf_i:psi_amp){
      //g :: g[sqrt][s_idx[mode_num]]
      //go through each spin
      complex<double> c_i = conf_i.second;
      double magnitude = abs(c_i);
      double pre_mul_mag = magnitude*dt;
      if(pre_mul_mag*g[max_level][s_idx[0]] < epsilon){
	break; //stop if c_i *H_max dt < epsilon
      }
      complex<double> pre_mul = -1i*c_i*dt;
      double reduced_epsilon = epsilon/pre_mul_mag;// dt g[][] c_i < ep -->g[][]<ep/(dt c_i)
      int j = 0;
      while((j < num_modes) && (g[max_level][s_idx[j]] >= reduced_epsilon) ){
	int level = conf_i.first.get_mode(j);
	
	/*for(int w = 0; w < psi_amp.size() ; w++){
	  for(int q = 0; q < num_modes; q++){
	    if( (psi_amp[w].first.get_mode(q))  > 0)
	      cout << psi_amp[w].first.get_mode(j) << endl;
	  }
	  }*/
	state_pair p;
	switch(level){
	case 0:
	  p.second = pre_mul*g[level+1][s_idx[j]];//only do the raising operator
	  p.first = conf_i.first;
	  p.first.increment_mode(j);
	  p.first.set_spin(~p.first.get_spin());
	  delta_psi.push_back(p);
	  break;
	case max_level:
	  p.second = pre_mul*g[level-1][s_idx[j]];//only do lowering op
	  p.first = conf_i.first;
	  p.first.decrement_mode(j);
  	  p.first.set_spin(~p.first.get_spin());
	  delta_psi.push_back(p);
	  break;
	default :
	  p.second = pre_mul*g[level+1][s_idx[j]];//do both
	  p.first = conf_i.first;
	  p.first.increment_mode(j);
	  p.first.set_spin(~p.first.get_spin());
	  delta_psi.push_back(p);

	  state_pair p2;
	  p2.second = pre_mul*g[level][s_idx[j]];
	  p2.first = conf_i.first;
	  p2.first.decrement_mode(j);
	  p2.first.set_spin(~p2.first.get_spin());
	  delta_psi.push_back(p);
	  break;
	}//end switch
	j++;
      }//end mode loop
    }//end config loop
  


    //~~~~~~~~~~merge and add existing configuration amplitudes~~~~~~~~~~~~~~~~~~~~~
    
    //delta psi is small, so do a small number of search ops to
    //do the duplicate problem beforehand
    typedef state_vec::iterator sIt;
    vector<pair<sIt,sIt>> eql_itors;
    //merge the duplicates between delta_psi and ps
    for(int i = 0; i<delta_psi.size(); i++ ){
      pair<sIt,sIt> eql= equal_range(psi_lbl.begin(),psi_lbl.end(),delta_psi[i],
				       [](auto &it1,auto&it2){return it1.first < it2.first ;}) ;
      if(eql.first !=eql.second){
	//found a dup
	eql_itors.push_back(eql);
	//add amplitude
	(eql.first)->second += delta_psi[i].second;
	//remove from delta_psi
	delta_psi.erase(delta_psi.begin()+i);
      }
    }
    
    //sort delta_psi
    sort(delta_psi.begin(),delta_psi.end(),
	 [](auto &it1,auto&it2){return it1.first < it2.first ;});
    //append delta_psi to end of psi 
    //reserve space so iterators don't change
    psi_lbl.reserve(psi_lbl.size() + delta_psi.size());
    //save midpoint for later
    state_vec::iterator midpoint_lbl = psi_lbl.end();
    //move the delta elements over to psi using move_iter, which
    //moves elements rather than copying them
    psi_lbl.insert(psi_lbl.end(),
		   make_move_iterator(delta_psi.begin()),
		   make_move_iterator(delta_psi.end()));
    
    inplace_merge(psi_lbl.begin(),midpoint_lbl,psi_lbl.end(),
		  [](auto &it1,auto&it2){return it1.first < it2.first ;});// [first,middle) and [middle,last)

    psi_amp.resize(psi_lbl.size());
    copy(psi_lbl.begin(),psi_lbl.end(),psi_amp.begin());
    sort(psi_amp.begin(),psi_amp.end(),[](auto it1,auto it2){
      return ( abs(it1.second) > abs(it2.second));});


    t +=dt;
    
  }

  

  return 0;
}
