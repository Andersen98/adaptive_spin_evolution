#include "hamiltonian.hpp"
#include "io_tools/input_tools.hpp"
hamiltonian::hamiltonian(const std::string &json_str):hamiltonian(load_json_str(json_str)){}

hamiltonian::hamiltonian(const param_vals &params_):g{{}},m{},params(params_),new2old{},old2new{},num_levels(1<<NUM_BITS),psi_delta(0),psi_lbl(0),psi_amp(0),state_connections(0),mode_cap_exceeded{}{

  
  vector<pair<int,double>> idx_g_pairs(NUM_MODES);
  for(int i=0;i< NUM_MODES;i++){
      idx_g_pairs[i] = std::make_pair(i,params.mode_couplings[i]);
    }
    
    sort(idx_g_pairs.begin(),idx_g_pairs.end(),[](auto &it1,auto &it2){return it1.second > it2.second;});

    //map from old idx to new idx
    //map new2old as well
    for(int i = 0; i < NUM_MODES; i++){
      int j =  idx_g_pairs[i].first;
      old2new[j] = i;
      new2old[i] = j;
    }

    for(int j = 0; j < NUM_MODES; j++){
	int i = idx_g_pairs[j].first;
	double gi = idx_g_pairs[j].second;
	//DEBUG
	assert(abs( gi - params.mode_couplings[i]) < .0001);
	//DEnd
	//descending order
	m[j] = params.mode_energies[i];
	
	for(int n = 0; n < num_levels; n++){
	  g[n][j] = gi*std::pow(n,.5);
	}
	//DEBUG
	assert(abs( g[1][j] - gi) < .0001);
	//DEnd
    }

    //construct initial states
    for(auto &sk: params.initial_state){
      int mode = old2new[sk.mode];
      state_ket k;
      k.set_mode(mode,sk.n);
      k.spin = sk.spin;
      k.amp = sk.amp;
      k.idx = state_ket::empty_idx;
      psi_lbl.push_back(k);
    }

    
    sort(psi_lbl.begin(),psi_lbl.end(),[](auto &it1,auto &it2){return it1<it2;});
    
    setup_connections();
    
    normalize_state(psi_lbl);
    psi_amp.resize(psi_lbl.size());
    copy(psi_lbl.begin(),psi_lbl.end(),psi_amp.begin());
    sort(psi_amp.begin(),psi_amp.end(),[](const state_ket &it1,const state_ket &it2){return norm(it1.amp) > norm(it2.amp);});
    
    
  }



  void hamiltonian::setup_connections(){

    //1 root node in psi_lbl
    psi_lbl[0].idx = state_connections.size();
    vector<dir_edge_mode> empty_edge(0);
    state_connections.push_back(empty_edge);
     
    for(auto &k: psi_lbl){
      if(k.idx == state_ket::empty_idx){
	vector<dir_edge_mode> edges;
	
	//apply hamiltonian
	for(int i = 0; i < NUM_MODES; i++){
	
	  ket_pair kp = get_connected_states(k,i);
	  bool raised = kp.raised.idx!=state_ket::null_idx;
	  bool lowered = kp.lowered.idx!=state_ket::null_idx;
	 
	if(raised){
	  //look for an instance
	  int vec_idx = binary_search_state(kp.raised,psi_lbl);
	  if(vec_idx>-1){
	    dir_edge_mode em;
	    em.out_idx = psi_lbl[vec_idx].idx;
	    em.connection_mode = i;
	    em.raised = true;
	    edges.push_back(em);
	  }
	  
	}
	if(lowered){
	  int vec_idx = binary_search_state(kp.lowered,psi_lbl);
	  if(vec_idx>-1){
	    dir_edge_mode em;
	    em.out_idx = psi_lbl[vec_idx].idx;
	    em.connection_mode = i;
	    em.raised = false;
	    edges.push_back(em); 
	  }
	}
	
      }//end mode loop

	k.idx = int(state_connections.size());
	state_connections.push_back(edges);
      
      }//end k loop

    }
    
  }

