#include "hamiltonian.hpp"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
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
  
  
  bool stop = false;
  for(int i = 0; i < int( psi_amp.size()); i++){
    stop = grow_configuration_space(i);
    if(stop){
      break;
	
    }
  }
  
  merge_states();
  
  append_connections(psi_amp);

  psi_lbl.resize(psi_amp.size());
  copy(psi_amp.begin(),psi_amp.end(),psi_lbl.begin());
  sort(psi_amp.begin(),psi_amp.end(),
       [](auto &it1,auto &it2){return norm(it1.amp) > norm(it2.amp);});


}

void hamiltonian::run_evolve(double dt){
  
  evolve_space(dt);
  normalize_state(psi_amp);
  
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

void  hamiltonian::run_step(complex<double> factor){
  //enter+leave with psi_amp set equal to psi_lbl set
  //enter+leave with psi_amp sorted by amplitude
  using namespace std;
  using namespace boost::numeric::ublas;
  v.resize(connection_matrix.size2());
  u.resize(connection_matrix.size2());
  for(uint i = 0; i < psi_lbl.size();i++){
    int idx = psi_lbl[i].idx;
    v[idx] = factor*psi_lbl[i].amp;
  }

  
  noalias(u) = v + factor*prod(connection_matrix,v);

  for(uint i = 0; i <psi_lbl.size();i++){
    int idx = psi_lbl[i].idx;
    psi_lbl[i].amp = u[idx];

  }
  
  normalize_state(psi_lbl);
  copy(psi_lbl.begin(),psi_lbl.end(),psi_amp.begin());
  sort(psi_amp.begin(),psi_amp.end(),[](const auto &it1,const auto &it2)
  {return std::norm(it1.amp)>std::norm(it2.amp);});
}

void hamiltonian::switch_evolve(){
  spin_up_matrix.resize(psi_lbl.size(),psi_lbl.size(),0,0);
  spin_down_matrix.resize(psi_lbl.size(),psi_lbl.size(),0,0);
  u.resize(connection_matrix.size2());
  v.resize(connection_matrix.size2());
  for(uint i=0; i < psi_lbl.size(); i++){
    int idx = psi_lbl[i].idx;
    v[idx] = psi_lbl[i].amp;
    if(psi_lbl[i].spin){
      spin_up_matrix(idx,idx) = 1;
      spin_down_matrix(idx,idx) = 0;
    }else{
      spin_up_matrix(idx,idx) = 0;
      spin_down_matrix(idx,idx) = 1;
    }
      
  }
  

  evolve_state = blas_state0;
  
  
}
void hamiltonian::blas_evolve(complex<double> factor){
  using namespace std;
  using namespace boost::numeric::ublas;

  switch(hamiltonian::evolve_state){
  case(blas_state0):
    noalias(u) = v + factor*prod(connection_matrix,v);
    u = (1.0/norm_2(u))*u;
    evolve_state = blas_state1;
    break;
  case(hamiltonian::blas_state1):
    noalias(v) = u + factor*prod(connection_matrix,u);
    v = (1.0/norm_2(v))*v;
    evolve_state = blas_state0;
    break;
  default:
    cout << "Error state is in grow, not blas_evolve. run switch_evolve()" <<endl;
    break;
  }
    
  
}

std::pair<double,double> hamiltonian::get_blas_spin_pop()const{
  
  using namespace std;
  using namespace boost::numeric::ublas;
  pair<double,double> result;


  switch(hamiltonian::evolve_state){
  case(blas_state0):
    result.first = norm_2_square(blas_vec(prod(spin_up_matrix,v)));
    result.second = norm_2_square(blas_vec(prod(spin_down_matrix,v)));
    break;
  case(hamiltonian::blas_state1):
    result.first = norm_2_square(blas_vec(prod(spin_up_matrix,u)));
    result.second = norm_2_square(blas_vec(prod(spin_down_matrix,u)));
    break;
  default:
    cout << "Error state is in grow, not blas_evolve. run switch_evolve()" <<endl;
    break;
  }
     
  return result;
}
