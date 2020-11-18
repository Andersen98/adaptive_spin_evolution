#include <vector>
#include <iostream>
#include <bitset>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <cmath>
#include <complex>
#include <utility>
#include <string>
#include <cassert>

#include "configuration.hpp"
#include "input_tools.hpp"
#include "hamiltonian.hpp"
#include "spin_density_matrix.hpp"

using namespace std;




int main(int argc, char * argv[]){
  
  //parameters to store system information
  //parameters also store simulation parameters like
  //energy tolerance and time step duration
  param_vals params;

  try{
    bool exit = get_params(params,argc,argv);
    if(exit)
      return 0;
  }catch(std::exception& e){
    cout << e.what() << endl;
    return 1;
  }
    

  
  const double final_time = params.tf;
  const double dt = params.dt;
  const int time_steps = final_time/dt;
  const double energy_cutoff = params.energy_cutoff;

 
  Spin_Params spin_params;
  spin_params.energy = params.atom_levels[1];

  //Complex Number Tools
  //also typedefs
  typedef complex<double> Amplitude;
  typedef bool spin_type;
  typedef State_Ket<spin_type,Amplitude, NUM_MODES,NUM_BITS> State_Ket;
  typedef vector<State_Ket> State_Vector;
  typedef vector<double>::iterator Iter_m; //mode iterator
  typedef vector<double>::iterator Iter_c; //coupling iterator
  //setup initial conditions
  State_Vector initial_state(1);
  initial_state[0].amp = 1;
  initial_state[0].spin = true;
  //make callback to calculate the population
  typedef Spin_Density_Matrix_Evolution<double> PC;
  PC matrix_recorder(dt,time_steps);
  //make hamiltonian
  hamiltonian<Iter_m,Iter_c,Spin_Params, State_Vector> h;
  
  h.setup(params.mode_energies.begin(),params.mode_couplings.begin(),spin_params,
      initial_state,energy_cutoff);



  //setup output stream
  string info_file = params.output_directory+to_string(params.run_id)+".info";
  string output_file = params.output_directory+to_string(params.run_id)+".out";
  string stats_file = params.output_directory+to_string(params.run_id)+".stats";
  std::ofstream of(output_file.c_str());
  std::ofstream of_stats(stats_file.c_str());
  std::ofstream of_info(info_file.c_str());
  if(of&& of_stats&&of_info){
    params.save(of_info);
    params.write_header(of);
    params.write_stats_header(of_stats);
    for(int i = 0; i <params.N; i++){
      State_Vector &v = h.get_psi_lbl();
      std::array<double,2> pop = matrix_recorder.get_spin_pop(v.begin(),v.end());
      params.write_pop_run(of,i+1,i*dt,pop[0],pop[1]);
      //    params.write_stats(of_stats, i+1, i*dt,v.size(), h.exceeded);
      h.simple_run(dt);
      cout << "run: (" << i << "/" << params.N<<") | ";
      cout << "number of states: " << v.size();
      cout.flush();
    }
  }
  cout <<endl;
  of.close();
  of_info.close();
  of_stats.close();
  
  return 0;
}
