
#include "hamiltonian.hpp"
#include <memory>
#include <type_traits>
using namespace std;

void hamiltonian::merge_states(){
    if(psi_delta.size()>0){
      state_vector::iterator psi_el;

      //sort cantidate states
      sort(psi_delta.begin(),psi_delta.end(),[](auto &it1, auto &it2){return it1 < it2; });
      
      //make space is psi_lbl for cantidate states
      psi_lbl.reserve(psi_amp.size()+psi_delta.size());


      //edge case: first n elements of psi_lbl are greater than
      
      
      psi_el = psi_amp.end();
      psi_amp.resize(psi_amp.size()+psi_delta.size());
     
      
     
      
    
      
      //dedupe
      vector<int> a(psi_delta.size());
      int j = 1;
      a[0] = 0;
      for(int i =1; i < psi_delta.size(); i++){
	if(psi_delta[a[j-1]] < psi_delta[i]){
	  a[j] = i;
	  j++;
	}
      }
      //copy over deduped list
      vector<state_ket> delta_clean(j);
      for(int i = 0; i < j; i++){

	delta_clean[i] = psi_delta[a[i]];
      }
      
      std::for_each(a.begin(),a.end(),[&](auto idx)
      state_ket w;
      vector<state_ket> delta_clean(psi_delta.size(),w);
      delta_clean[0] = psi_delta[0];
      delta_clean[0].amp = complex<double>(0,0);
      int j = 1;
      for(int i = 1; i < psi_delta.size(); i++){
	if(delta_clean[j-1] < psi_delta[i]){
	  delta_clean[j] = psi_delta[i];
	  delta_clean[j].amp = complex<double>(0,0);
	  j++;

	}	
      }
      //RESIZE DELTA CLEAN BEFORE DOING ANYTHING ELSE!
      int delta_clean_size = j;
      delta_clean.resize(delta_clean_size,w);      

      psi_el = std::set_difference(delta_clean.begin(),delta_clean.end(),
				   psi_lbl.begin(),psi_lbl.end(),psi_delta.begin(),[](auto &it1, auto &it2){return it1 < it2; });

      
      psi_delta.resize(std::distance(psi_delta.begin(),psi_el));
      //append new connections before merge
      append_connections(psi_delta);


      //merge into psi_lbl
      psi_amp.resize(delta_clean.size()+psi_lbl.size(),w);    
      psi_el = std::set_union(psi_lbl.begin(),psi_lbl.end(),psi_delta.begin(),psi_delta.end(),
			      psi_amp.begin(),[](auto &it1, auto &it2){return it1 < it2; });


      psi_amp.resize(std::distance(psi_amp.begin(),psi_el));


      
    }//endif delta_size.()
    
  
    
  }
