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
void hamiltonian::append_connections(state_vector &psi_mixed ){


  if(psi_mixed.size() > connection_matrix.size1()){
    grow_matrix(2*psi_mixed.size());
    
  }

  for(auto &ket : psi_mixed){
    
    if(ket.idx == state_ket::empty_idx){
      
      uint update_idx = state_connections.size();
      update_matrix_diag_only_spin(update_idx, ket.spin);
      vector<dir_edge_mode> edges;   
      //apply hamiltonian
      for(int i = 0; i < NUM_MODES; i++){
	 
	ket_pair kp = get_connected_states(ket,i);
	bool raised = kp.raised.idx==state_ket::empty_idx;
	bool lowered = kp.lowered.idx==state_ket::empty_idx;
	uint level = ket.get_mode(i);
	update_matrix_diag_only_mode(update_idx,i,level);
	if(raised){
	  //look for an instance
	  int vec_idx = binary_search_state(kp.raised,psi_mixed);
	  if(vec_idx > -1){
	    dir_edge_mode em;
	    em.out_idx = psi_mixed[vec_idx].idx;
	    em.connection_mode = i;
	    em.raised = true;
	    edges.push_back(em);
	    connection_matrix(update_idx,em.out_idx) =
	      g[level+1][i];
	  }
	}
	if(lowered){
	  int vec_idx = binary_search_state(kp.lowered,psi_mixed);
	  if(vec_idx>-1){
	    dir_edge_mode em;
	    em.out_idx = psi_mixed[vec_idx].idx;
	    em.connection_mode = i;
	    em.raised = false;
	    edges.push_back(em);
	    connection_matrix(update_idx,em.out_idx) =
	      g[level][i];
	  }
	}

      }//end mode loop

    
      //'color' k
      ket.idx = int(state_connections.size());

      
      state_connections.push_back(edges);
    }//end if empty_idx
  }//end k loop
}
