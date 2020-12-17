
#include "hamiltonian.hpp"


void hamiltonian::evolve_space(double dt){



    vector<int> con2amp(psi_amp.size());
    for(int i = 0; i < psi_amp.size(); i++){
      //con -> ampidx
      con2amp[psi_amp[i].idx] = i;
    }
    vector<complex<double>> delta(psi_amp.size());

    for(int i = 0; i < state_connections.size(); i ++){
      complex<double> start_amp = psi_amp[con2amp[i]].amp;
      //diagonal term
      //--atom energy
      delta[i] += start_amp*complex<double>(0,-1)*dt*(psi_amp[con2amp[i]].spin? params.up_energy:params.down_energy);

      assert(i==psi_amp[con2amp[i]].idx);
      for(auto & edge:state_connections[i]){
	int out_idx = edge.out_idx;
	complex<double> out_amp = psi_amp[con2amp[out_idx]].amp;
	int connection_mode = edge.connection_mode;
	int start_level = psi_amp[con2amp[i]].get_mode(connection_mode);
	bool raised =  edge.raised;
	//start -> finish
	double g_factor = raised?(g[start_level+1][connection_mode]):(g[start_level][connection_mode]);
	delta[out_idx] += start_amp*complex<double>(0,-1)*dt*g_factor;
	//diagonal term
	delta[i] += start_amp*complex<double>(0,-1)*dt*double(start_level)*m[connection_mode];

	//finish->start
	delta[i] += out_amp*complex<double>(0,-1)*dt*g_factor;

      }
    }


    for(int i = 0; i < psi_amp.size(); i++){
      psi_amp[con2amp[i]].amp += delta[i];
      
    }



}
