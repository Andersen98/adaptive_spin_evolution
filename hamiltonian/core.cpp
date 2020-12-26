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
    N += std::norm(k.amp);
  }
  N = 1/std::sqrt(N);
  std::transform(p.begin(),p.end(),p.begin(),
		 [&](auto k){ k.amp = N*k.amp;return k;}); 
  
}

void hamiltonian::update_matrix_diag(const state_ket &k){
  

  update_matrix_diag_only_spin((uint)k.idx,k.spin);
  
  for(int i = 0; i < NUM_MODES; i++){
    
    int level = k.get_mode(i);
    update_matrix_diag_only_mode((uint)k.idx,i,level);
  }

}

void hamiltonian::update_matrix_diag_only_spin(uint update_idx, bool spin){

  switch(spin){
  case(true):
    connection_matrix(update_idx,update_idx) += params.up_energy;
    break;
  default:
    connection_matrix(update_idx,update_idx) += params.down_energy;
    break;
    
  }
}
void hamiltonian::update_matrix_diag_only_mode(uint update_idx, int mode, int level){

  connection_matrix(update_idx,update_idx) += (double)level*m[mode];
  
  
}
void hamiltonian::grow_matrix(uint new_size){
  assert(new_size > connection_matrix.size1());
  base_matrix_type b( boost::numeric::ublas::zero_matrix<double>(new_size,new_size));
  matrix_type c(b);
  for(uint i=0;i < connection_matrix.size1(); i++){
    for(uint j=0; j <=i; j++){
      c(i,j) = connection_matrix(i,j);
    }
  }
  base_matrix = b;
  connection_matrix = c;

}
