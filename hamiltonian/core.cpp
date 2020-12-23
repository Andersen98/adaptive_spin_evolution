#include "hamiltonian.hpp"


hamiltonian::ket_pair::ket_pair():raised(),lowered(){}
hamiltonian::ket_pair::ket_pair(const state_ket &k){
  raised = k;
  lowered = k;
  raised.amp = 0;
  lowered.amp = 0;
}


hamiltonian::ket_pair hamiltonian::get_connected_states(const state_ket &k,const int mode){

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
     case((1<<NUM_BITS)-1):
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




void hamiltonian::normalize_state(state_vector &p){
  
  double N = 0;
  for(const state_ket &k: p){
    N += norm(k.amp);
  }
  N = 1/std::sqrt(N);
  for_each(p.begin(),p.end(),[N](state_ket &k){k.amp*=N;});
  
}
