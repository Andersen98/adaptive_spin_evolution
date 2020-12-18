#include "hamiltonian.hpp"


using namespace std;

const hamiltonian::state_vector &hamiltonian::get_state_vector()const{
    return psi_lbl;
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

pair<double,double> hamiltonian::get_emitter_cavity_prob(bool emitter, int cav_n)const{
  pair<double,double> result;
  pair<state_ket,state_ket> kp;
  state_ket k1,k2;
  k1.set_mode(0,cav_n);
  k1.spin = emitter;

  if(emitter){
    //k1 is emmitter
    //k2 is cav
    //corresponding stae is spin down +1 cav
    k2.set_mode(0,cav_n+1);
    k2.spin = false;
    
    kp.first = k1;
    kp.second = k2;
  }else{
    //k2 is emmiter
    //k1 is cav
    assert( cav_n > 0);
    k2.set_mode(0, cav_n-1);
    k2.spin = true;

    kp.first = k2;
    kp.second = k1;
  }

  int idx1,idx2;
  idx1 = binary_search_state(kp.first,psi_lbl);
  idx2 = binary_search_state(kp.second,psi_lbl);
  for(const state_ket &k:psi_lbl){
    cout << k << endl;
  }
  cout <<idx1 << '\t' << idx2 << endl;
  result.first = norm(psi_lbl[idx1].amp);
  result.second = norm(psi_lbl[idx2].amp);

  cout << psi_lbl[idx1] << endl;
  cout << psi_lbl[idx2] << endl;
  

  return(result); 


}

array<tuple<int,double,double>,NUM_MODES> hamiltonian::get_modeLbl_quanta_pop()const{
    array<tuple<int,double,double>,NUM_MODES> result;

    for_each(result.begin(),result.end(),[j=0](auto &el)mutable{el= std::make_tuple<int,double,double>(j++,0,0);});
    for_each(result.begin(),result.end(),[&](auto &el){
      for(auto &ket:psi_lbl){
	double probability =norm(ket.amp);
	int number_expectation =0;
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
  
