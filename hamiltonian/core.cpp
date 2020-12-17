#include "hamiltonian.cpp"

hamiltonian::void run_step(double dt){
    mode_cap_exceeded.clear();
    psi_delta.clear();
    psi_delta.reserve(NUM_MODES*psi_amp.size());

    //TODO: Maybe add a threshold for number of runs where the
    //configuration has not grown. 
    bool stop = false;
    for(int i = 0; i < psi_amp.size(); i++){
      stop = grow_configuration_space(i);
      if(stop){
	break;
	
      }
    }
    
    merge_states();

    evolve_current_space(dt);
    
  }

ket_pair hamiltonian::get_connected_states(const state_ket &k,const int mode){

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



static double hamiltonian::abs_sqrd(const state_ket &p){
  return(std::norm(p.amp));
}

static double hamiltonian::abs(const state_ket &p){
  return(std::abs(p.amp));
}

static void hamiltonian::normalize_state(state_vector &p){
  
  double N = 0;
  for(const state_ket &k: p){
    N += abs_sqrd(k);
  }
  N = 1/std::sqrt(N);
  for_each(p.begin(),p.end(),[N](state_ket &k){k.amp*=N;});
  
}
