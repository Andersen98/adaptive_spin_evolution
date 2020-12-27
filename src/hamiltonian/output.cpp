#include "hamiltonian.hpp"


using namespace std;

hamiltonian::state_vector hamiltonian::get_state()const{

  return state_vector(psi_lbl);
}

hamiltonian::matrix_type &hamiltonian::get_matrix(){
  return(connection_matrix);
}
int hamiltonian::get_psi_size()const{
    return psi_lbl.size();
  }
    
  //returns current spin up, spin down 
pair<double,double> hamiltonian::get_spin_pop()const{
    pair<double,double> result = std::make_pair(0,0);

    for(const auto&k:psi_lbl){
      if(k.spin){
	result.first += norm(k.amp);
      }else{
	result.second += norm(k.amp);
      }	
    }
    
    return result;

}

pair<double,double> hamiltonian::get_emitter_cavity_prob()const{
  pair<double,double> result;
  int num_quanta = -1;
  //look for initial state (idx 0 and 1 are assumed to be on same N manifold )
  for(const auto &k:psi_lbl){

    if(k.idx >1){
      continue;
    }else{
      assert( (k.idx==1||k.idx==0));
      int mode_level = k.get_mode(0);
      if(k.spin){
	result.first = norm(k.amp);
	if(num_quanta >-1){
	  assert( (mode_level+1 == num_quanta));
	}
	num_quanta = 1 + mode_level;
      }else{
	result.second = norm(k.amp);
	if(num_quanta > -1){
	  assert( mode_level == num_quanta);
	  
	}
      }

    }
    
  }
  
  return(result); 


}

array<tuple<int,double,double>,NUM_MODES> hamiltonian::get_modeLbl_quanta_pop()const{
    array<tuple<int,double,double>,NUM_MODES> result;

    for_each(result.begin(),result.end(),[j=0](auto &el)mutable{el= std::make_tuple<int,double,double>(j++,0,0);});
    for_each(result.begin(),result.end(),[&](auto &el){
      for(auto &ket:psi_lbl){
	double probability =norm(ket.amp);
	int modeLbl = std::get<0>(el);
	int nVal = ket.get_mode(modeLbl);
	std::get<1>(el) += nVal*probability;
	std::get<2>(el) += probability;
	  
      };});
    return result;
  }

std::array<int,NUM_MODES> hamiltonian::get_new2old()const{
  array<int,NUM_MODES> result = {};
  copy(new2old.begin(),new2old.end(),result.begin());
  return(result);
  

}
std::array<int,NUM_MODES> hamiltonian::get_old2new()const{
  array<int,NUM_MODES> result = {};
  copy(old2new.begin(),old2new.end(),result.begin());
  return(result);

}
  
