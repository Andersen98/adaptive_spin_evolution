#include "hamiltonian.hpp"

int hamiltonian::binary_search_lbl(const state_ket &k){
    std::pair<state_vector_iterator,state_vector_iterator> psi_el;
    int result = -1; //if it fails, the result is -1;
    psi_el = std::equal_range(psi_lbl.begin(),psi_lbl.end(),k,
			      [](auto &el,auto &val){return el < val;});

    if(psi_el.first != psi_el.second){
      result = std::distance(psi_lbl.begin(),psi_el.first);
    }

    return(result);
  }

void hamiltonian::append_connections(){

   
    
  for(auto &ket : psi_amp){
    
    if(ket.idx != state_ket::empty_idx){
      continue;
    }
    
    ket.idx = int(state_connections.size());
    vector<dir_edge_mode> edges;   
    //apply hamiltonian
    for(int i = 0; i < NUM_MODES; i++){
	 
      ket_pair kp = get_connected_states(ket,i);
      bool raised = kp.raised.idx==state_ket::empty_idx;
      bool lowered = kp.lowered.idx==state_ket::empty_idx;

      if(raised){
	//look for an instance
	int vec_idx = binary_search_lbl(kp.raised);
	if(vec_idx > -1){
	  dir_edge_mode em;
	  em.out_idx = psi_amp[vec_idx].idx;
	  em.connection_mode = i;
	  em.raised = true;
	  edges.push_back(em);
	}
      }
      if(lowered){
	int vec_idx = binary_search_lbl(kp.lowered);
	if(vec_idx>-1){
	  dir_edge_mode em;
	  em.out_idx = psi_amp[vec_idx].idx;
	  em.connection_mode = i;
	  em.raised = false;
	  edges.push_back(em); 
	}
      }

    }//end mode loop
    
    state_connections.push_back(edges);
       
  }//end k loop
}
     
