#include "hamiltonian.hpp"
#include <unsupported/Eigen/MatrixFunctions>
void hamiltonian::set_epsilon(double e){
  
  params.energy_cutoff = e;

}

void hamiltonian::run_grow(){
  //assumes psi_amp is sorted by amplitude entering and leaving
  //assumes psi_lbl is sorted by label  entering and leaving
  
  mode_cap_exceeded.fill(-1);
  psi_delta.resize(0);
  psi_delta.reserve(NUM_MODES*psi_amp.size());
  
  
  bool stop = false;
  for(int i = 0; i < int( psi_amp.size()); i++){
    stop = grow_configuration_space(i);
    if(stop){
      break;
	
    }
  }
  
  merge_states();
  
  for(auto &k: psi_amp){
    add_connection(k,psi_amp);
  }
 
  psi_lbl.resize(psi_amp.size());
  copy(psi_amp.begin(),psi_amp.end(),psi_lbl.begin());
  sort(psi_amp.begin(),psi_amp.end(),
       [](auto &it1,auto &it2){return norm(it1.amp) > norm(it2.amp);});


}


void hamiltonian::set_zero_except_init(){
  using namespace std;
  int init_size = params.initial_state.size();

  double init_norm = 0;
  for(const auto &k_init:params.initial_state){
    init_norm += std::norm(k_init.amp);
  }
  init_norm = std::sqrt(init_norm);
  
  for( auto &k:psi_lbl){

    if(k.idx < init_size){
      //this part of initial state
      for(auto &ik:params.initial_state){
	if(k == ik){
	  k.amp =(1.0/init_norm) *ik.amp;
	}
      }
    }else{
      k.amp = 0;
    }
  }

  copy(psi_lbl.begin(),psi_lbl.end(),psi_amp.begin());
  sort(psi_amp.begin(),psi_amp.end(),[](const auto &it1,const auto &it2){
    return norm(it1.amp)>norm(it2.amp);
  });
    
}
void hamiltonian::store_vector(){
  psi_u.resize(psi_lbl.size());
  psi_uinit.resize(psi_lbl.size());
  SpinMatrix.resize(2,psi_lbl.size());
  for(uint i = 0; i < psi_lbl.size(); i++){
    int idx = psi_lbl[i].idx;
    psi_u[idx] = psi_lbl[i].amp;
    if(idx < params.initial_state.size()){
      int j = binary_search_state(psi_lbl[i],params.initial_state);
      psi_uinit[idx] = params.initial_state[j].amp;
    }else{
      psi_uinit[idx] = 0;
    }

    if(psi_lbl[i].spin){
      SpinMatrix(0,idx) =1;
      SpinMatrix(1,idx) = 0;
    }else{
      SpinMatrix(0,idx) = 0;
      SpinMatrix(1,idx) = 1;
    }
  }
 
  
}
void hamiltonian::store_matrix(){
  using namespace Eigen;
  H_matrix.resize(psi_lbl.size(),psi_lbl.size());
  H_matrix.setFromTriplets(state_connections.begin(),state_connections.end());
  H_exp.resize(psi_lbl.size(),psi_lbl.size());
  
}
std::pair<double,double> hamiltonian::evolve_state(complex<double> time){
  std::pair<double,double> result;
  SpMat U = H_exp.pow(complex<double>(0,-time));
  Eigen::Matrix<std::complex<double>, 2, Eigen::Dynamic> v = SpinMatrix*(U*psi_uinit);
  Eigen::Matrix<std::complex<double>, 2, 1> u = v.colwise().norm();
  result = std::make_pair(u(0,0),u(0,1));

}
void  hamiltonian::run_step(complex<double> factor){
  //enter+leave with psi_amp set equal to psi_lbl set
  //enter+leave with psi_amp sorted by amplitude
  SpMat A(next_idx,next_idx);
  A.setFromTriplets(state_connections.begin(),state_connections.end());
  ComplexVec b(psi_lbl.size());

  for(uint i = 0; i < psi_lbl.size(); i++){
    int idx = psi_lbl[i].idx;
    b[idx] = psi_lbl[i].amp;
  }
  //did not expect such an awkward expression (sorry)
  b = b + A.selfadjointView<Eigen::Upper>()*(factor*b);

  for(uint i = 0; i < psi_lbl.size(); i++){

    int idx = psi_lbl[i].idx;
    psi_lbl[i].amp = b[idx];
    
  }
  normalize_state(psi_lbl);
  std::copy(psi_lbl.begin(),psi_lbl.end(),psi_amp.begin());
  std::sort(psi_amp.begin(),psi_amp.end(),
	    [](const auto &it1,const auto &it2){return std::norm(it1.amp) > std::norm(it2.amp);});
  
}

void hamiltonian::reset_with_state(state_vector v){
  std::copy(v.begin(),v.end(),params.initial_state.begin());
  reset();
}
void hamiltonian::reset(){

  //construct initial states
  psi_delta.resize(0);
  psi_lbl.resize(0);
  psi_amp.resize(0);
  state_connections.resize(0);
  next_idx = 0;
  for(auto &k: params.initial_state){
    state_ket k_transform;
      k_transform.idx = state_ket::empty_idx;
      k_transform.amp = k.amp;
      k_transform.spin = k.spin;
      for(uint i = 0; i < NUM_MODES; i++){
	int level = k.get_mode(new2old[i]);
	k_transform.set_mode(i,level);
      }
      psi_lbl.push_back(k_transform);
    }

    
    sort(psi_lbl.begin(),psi_lbl.end(),[](auto &it1,auto &it2){return it1<it2;});
    
    setup_connections();
    
    normalize_state(psi_lbl);
    psi_amp.resize(psi_lbl.size());
    copy(psi_lbl.begin(),psi_lbl.end(),psi_amp.begin());
    sort(psi_amp.begin(),psi_amp.end(),[](const state_ket &it1,const state_ket &it2){return norm(it1.amp) > norm(it2.amp);});
    


}
