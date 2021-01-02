#include "hamiltonian.hpp"

//search sorted vector for a given ket k
int hamiltonian::binary_search_state(const state_ket &k,const state_vector &psi_search)const{

    int result = -1; //if it fails, the result is -1;
    auto psi_el = std::equal_range(psi_search.begin(),psi_search.end(),k,
			      [](const state_ket &el,const state_ket &val){return el < val;});

    if((psi_el.first != psi_el.second)&& (psi_el.first->idx != state_ket::empty_idx )  ){
      result = std::distance(psi_search.begin(),psi_el.first);
    }

    return(result);
  }



//psi_mixed is a sorted vector where some elements are 'colored' empty
//and other elements are 'colored' with a defined idx.
//each step, we find an empy element 'k', find k's connections
//to non empy elements. then 'color' k with an idx
//-----------------------
//adds a column to sparse matrix for a given state @ket
//here row < column so this constructs an upper triangular
//matrix
void hamiltonian::add_connection(state_ket &ket,state_vector &psi_mixed){
  if(ket.idx == state_ket::empty_idx){
    complex<double> diag_val = 0;
    diag_val += matrix_diag_only_spin(ket.spin);
    int update_idx = next_idx++;
    
  
    //apply hamiltonian
    for(int i = 0; i < NUM_MODES; i++){	 
      ket_pair kp = get_connected_states(ket,i);
      bool raised = kp.raised.idx==state_ket::empty_idx;
      bool lowered = kp.lowered.idx==state_ket::empty_idx;
      int level = ket.get_mode(i);
      diag_val += matrix_diag_only_mode(i,level);
      if(raised){
	//look for an instance
	int vec_idx = binary_search_state(kp.raised,psi_mixed);
	if(vec_idx > -1){
	  int row = psi_mixed[vec_idx].idx; //desination
	  int col = update_idx;
	  complex<double> val = complex<double>(g[level+1][i],0);
	  state_connections.push_back(T(row,col,val));
	}
      }
      if(lowered){
	int vec_idx = binary_search_state(kp.lowered,psi_mixed);
	if(vec_idx>-1){
	  int row = psi_mixed[vec_idx].idx;
	  int col = update_idx;
	  complex<double> val = complex<double>(g[level][i],0);
	  state_connections.push_back(T(row,col,val));
	}
      }
      
    }//end mode loop
    state_connections.push_back(T(update_idx,update_idx,diag_val));
    //'color' k
    ket.idx = update_idx; 
    
  }//end if
}
  

