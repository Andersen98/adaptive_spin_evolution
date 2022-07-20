#include "hamiltonian.hpp"
#include "../io_tools/input_tools.hpp"


hamiltonian::hamiltonian(const std::string &json_str):hamiltonian(load_json_str(json_str)){}

hamiltonian::hamiltonian(const param_vals &params_):g{{}},m{},params(params_),new2old{},old2new{},num_levels(1<<NUM_BITS),psi_delta(0),psi_lbl(0),psi_amp(0),state_connections(0),next_idx(0),mode_cap_exceeded{}{

  
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
      g[n][j] = gi*std::sqrt(n);
    }
    //DEBUG
    assert(abs( g[1][j] - gi) < .0001);
    //DEnd
  }
  reset();
}



void hamiltonian::setup_connections(){
    
  //1 root node in psi_lbl
  psi_lbl[0].idx = next_idx++;
  std::complex<double> diag_val = matrix_diag(psi_lbl[0]);
  state_connections.push_back(T(psi_lbl[0].idx,psi_lbl[0].idx,diag_val));
  
  for(auto &k: psi_lbl){
      add_connection(k, psi_lbl);
  }				
}


