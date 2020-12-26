#include "hamiltonian.hpp"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

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
  typedef boost::numeric::ublas::vector<complex<double>> blas_vec;
  blas_vec v(connection_matrix.size2());
  for(uint i = 0; i < psi_lbl.size();i++){
    int idx = psi_lbl[i].idx;
    v[idx] = factor*psi_lbl[i].amp;
  }

  
  blas_vec u(boost::numeric::ublas::prod(v,connection_matrix));
  
  for(uint i = 0; i <psi_lbl.size();i++){
    int idx = psi_lbl[i].idx;
    psi_lbl[i].amp += u[idx];
  }
  
  normalize_state(psi_lbl);
  copy(psi_lbl.begin(),psi_lbl.end(),psi_amp.begin());
  sort(psi_amp.begin(),psi_amp.end(),[](const auto &it1,const auto &it2)
  {return std::norm(it1.amp)>std::norm(it2.amp);});
}
