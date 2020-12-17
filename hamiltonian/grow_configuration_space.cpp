#include "hamiltonian.hpp"


//returns true to halt evaluation
bool hamiltonian::grow_configuration_space(int idx){
    bool stop = false;

    double magnitude = abs(psi_amp[idx].amp);
    if(g[num_levels-1][0]*magnitude < params.energy_cutoff){
      return true;
    }

    //diagonal term
    // Sum omega_j n_j + params.emitter_energy;
    int j = 0;
    while((j < NUM_MODES) && (params.energy_cutoff < magnitude*g[num_levels-1][j]) ){

      ket_pair kp = get_connected_states(psi_amp[idx],j);
      int level = psi_amp[idx].get_mode(j);

      
      if( (kp.raised.idx != state_ket::null_idx) && (magnitude* g[level+1][j]> params.energy_cutoff) ){
	psi_delta.push_back(kp.raised);
      }
      if( (kp.lowered.idx != state_ket::null_idx) && (magnitude*g[level+1][j] > params.energy_cutoff)){
	psi_delta.push_back(kp.lowered);
      }	
      j++;
    }
   
    return(stop);

}
