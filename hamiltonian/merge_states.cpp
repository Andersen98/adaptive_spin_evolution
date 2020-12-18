#include "hamiltonian.hpp"
#include <memory>
#include <type_traits>
#include <iterator>
using namespace std;

void hamiltonian::merge_states(){
    if(psi_delta.size()>0){
      state_vector::iterator psi_mid,psi_end;
      
      //sort cantidate states
      sort(psi_delta.begin(),psi_delta.end(),
 	   [](const state_ket &it1,const state_ket& it2){return it1 < it2; });

      
      
      

      state_vector delta_clean(psi_delta.size());
      //append unique copy of psi_delta to end of psi_lbl
      state_vector::iterator delta_end =
	unique_copy(psi_delta.begin(),psi_delta.end(),
		    delta_clean.begin(),
		    [](const state_ket &it1,const state_ket &it2)
		    {return it1==it2;});
      delta_clean.resize(distance(delta_clean.begin(),delta_end));
      
      
      psi_amp.resize(psi_lbl.size()+delta_clean.size());
      //sort the two sequences
      auto amp_end = set_union( psi_lbl.begin(),psi_lbl.end(),delta_clean.begin(),delta_clean.end(),psi_amp.begin(),
		     [](const state_ket &it1,const state_ket& it2){return it1 < it2; });
      psi_amp.resize(distance(psi_amp.begin(),amp_end));
	    
      
    }//endif delta_size.()
    
  
    
  }
