#include "hamiltonian.hpp"
void par_test_two(){}
void hamiltonian::set_epsilon(double e){

  params.energy_cutoff = e;

}
void hamiltonian::run_grow_evolve(double dt){
  mode_cap_exceeded.fill(-1);
  psi_delta.resize(0);
  psi_delta.reserve(NUM_MODES*psi_amp.size());

  //TODO: Maybe add a threshold for number of runs where the
  //configuration has not grown. 
  bool stop = false;
  for(int i = 0; i < int( psi_amp.size()); i++){
    stop = grow_configuration_space(i);
    if(stop){
      break;
	
    }
  }
  
  merge_states();
  
  append_connections(psi_amp);
  
  evolve_space(dt);
  
  normalize_state(psi_amp);
  psi_lbl.resize(psi_amp.size());
  copy(psi_amp.begin(),psi_amp.end(),psi_lbl.begin());
  sort(psi_amp.begin(),psi_amp.end(),
       [](auto &it1,auto &it2){return norm(it1.amp) > norm(it2.amp);});
}
void hamiltonian::run_grow(){
  //assumes psi_amp is sorted by amplitude entering and leaving
  //assumes psi_lbl is sorted by label  entering and leaving
  mode_cap_exceeded.fill(-1);
  psi_delta.resize(0);
  psi_delta.reserve(NUM_MODES*psi_amp.size());
  
  //TODO: Maybe add a threshold for number of runs where the
  //configuration has not grown. 
  bool stop = false;
  for(int i = 0; i < int( psi_amp.size()); i++){
    stop = grow_configuration_space(i);
    if(stop){
      break;
	
    }
  }
  
  merge_states();
  
  append_connections(psi_amp);

  psi_lbl.reserve(psi_amp.size());
  copy(psi_amp.begin(),psi_amp.end(),psi_lbl.begin());
  sort(psi_amp.begin(),psi_amp.end(),
       [](auto &it1,auto &it2){return norm(it1.amp) > norm(it2.amp);});


}

void hamiltonian::run_evolve(double dt){
  
  evolve_space(dt);
  normalize_state(psi_amp);
  
}

void hamiltonian::set_zero_except_init(){

  int init_size = params.initial_state.size();

  double init_norm = 0;
  for(const auto &k_init:params.initial_state){
    init_norm += std::norm(k_init.amp);
  }
  init_norm = std::sqrt(init_norm);
  
  for( auto &k:psi_lbl){

    if(k.idx < init_size){
      //this part of initial state
      k.amp =(1.0/init_norm) *params.initial_state[k.idx].amp;
    }else{
      k.amp = 0;
    }
  }

  
  
  
}
