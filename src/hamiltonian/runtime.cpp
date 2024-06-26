#include "hamiltonian.hpp"
#include <omp.h>
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
/* TODO FIX THIS
 * vector<complex<double>> hamiltonian::get_decay_diags(int delta_e){

  vector<complex<double>> decay_diags(psi_lbl.size());
  
  //for each vector in psi, get the connected states
  for(int i=0; i < int(psi_amp.size()); i++){
    double E0 = std::real(matrix_diag(psi_amp[i]));
    int idx = psi_amp[i].idx;
    for(int j = 0; j < NUM_MODES; j++){
      ket_pair kp = get_connected_states(psi_amp[i],j);
      int level = psi_amp[i].get_mode(j);
      if( (kp.raised.i != state_ket::null_i) && (magnitude* g[level+1][j]> params.energy_cutoff) ){
	//Check to see if energy is within delta_e
	double E_k = std::real(matrix_diag(kp.raised ));
	double V_ik = g[level+1][j];
	if(std::abs(E_k-E0) <= delta_e){
	  if(binary_search_state(kp.raised, psi_lbl) < 0){
	    //state is not in current space, so we add it as decay term
	    decay_diags[idx] += complex<double>(0,-(V_ik*V_ik)/delta_e);
	  }
	}
      }
      if( (kp.lowered.i != state_ket::null_i) && (magnitude*g[level][j] > params.energy_cutoff)){
	//Check to see if energy is within delta_e
	double E_k = std::real(matrix_diag(kp.lowered ));
	double V_ik = g[level][j];
	if(std::abs(E_k-E0) <= delta_e){
	  //state is within energy threshold. 
	  if(binary_search_state(kp.lowered,psi_lbl) < 0){
	    //State is not in current space, add it as complex energy
	    //Given by fermi's golden rule
	    decay_diags[idx] += complex<double>(0,(V_ik*V_ik)/delta_e);
	  }
	}
      }

    }
  }
  
  return(decay_diags);
}
*/

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
      state_ket k = psi_lbl[i];
      int j = binary_search_state(k,params.initial_state);
      assert(j >=0);
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

  

  
}
struct Exp{
  std::complex<double> operator()(std::complex<double> x)const{
    return std::exp(x);
  }    

};

std::pair<double,double> hamiltonian::evolve_state(double time){
  using namespace Eigen;
  setNbThreads(12); 
  std::pair<double,double> result;
  std::cout << nbThreads()<<std::endl;

  ComplexVec eigen_val_exp = (std::complex<double>(0,-time)*H_eigen_vals).unaryExpr(Exp());
  
  ComplexVec u = H_eigen_vectors*eigen_val_exp.asDiagonal()*H_eigen_vectors.adjoint()*psi_uinit;



  //  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> u = SpinMatrix*( complex<double>(0,-time)*H_exp).exp()*psi_uinit;
  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> v = u.array()*u.conjugate().array();
  ComplexVec u2 = SpinMatrix*v;
  result = std::make_pair(std::real(u2(0,0)),std::real(u2(1,0)));
  //  std::cout <<u2 << std::endl;
  return result;
}

std::vector<std::tuple<int,int,std::complex<double>>> hamiltonian::get_tuples()const{
  std::vector<std::tuple<int,int,std::complex<double>>>  result;
  using namespace std;
  for(auto &t:state_connections){

    result.push_back(tuple<int,int,complex<double>>{t.row(),t.col(),t.value()});
  }
  
  return result;
}
std::vector<std::pair<int, int>> hamiltonian::get_spin_idxs()const{
  using namespace std;
  vector<pair<int, int> > result(psi_lbl.size());
  for(const auto &k: psi_lbl){
    if(k.spin){
      result[k.idx] = make_pair(1,0);
    }else{

      result[k.idx] = make_pair(0,1);
    }
    
  }
  return result;

}

std::pair<double,double>  hamiltonian::evolve_step(complex<double> factor){
  
  using namespace std;
  using namespace Eigen;
  omp_set_num_threads(16);
  psi_u = psi_u + factor*(H_matrix.selfadjointView<Upper>()*psi_u).eval();
  psi_u.normalize();
  pair<double,double> result;
  ComplexVec u = psi_u.conjugate().array()*psi_u.array();
  u = SpinMatrix*u;
  

  result.first = real(u[0]);
  result.second = real(u[1]);
  return result;
  

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
  std::sort(params.initial_state.begin(),params.initial_state.end(),[](auto &it1, auto &it2){return it1 < it2;});
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
