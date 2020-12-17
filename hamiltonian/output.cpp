#include "hamiltonian.hpp"


using namespace std;

state_vector &hamiltonian::get_state_vector()const{
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
	result.first += abs_sqrd(k);
      }else{
	result.second += abs_sqrd(k);
      }	
    }
    
    return result;

  }

array<tuple<int,double,double>,NUM_MODES> hamiltonian::get_modeLbl_quanta_pop(){
    array<tuple<int,double,double>,NUM_MODES> result;

    for_each(result.begin(),result.end(),[j=0](auto &el)mutable{el= std::make_tuple<int,double,double>(j++,0,0);});
    for_each(result.begin(),result.end(),[&](auto &el){
      for(auto &ket:psi_lbl){
	double probability =abs_sqrd(ket);
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
  
