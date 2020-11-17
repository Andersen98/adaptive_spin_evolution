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
  //const double high_freq = 60*1.81;
  const double high_freq = params.largest_frequency;
  const double dt = 1/high_freq;
  const int time_steps = final_time/dt;
  params.N = time_steps;
  const double energy_cutoff = params.energy_cutoff;
  const int num_spins = 2;
  const int num_modes = 1000;
  const int num_bits = 8;
  const long int max_level = (1<<num_bits) -1;

 
  Spin_Params spin_params;
  spin_params.energy = params.atom_levels[1];

  //Complex Number Tools
  //also typedefs
  typedef complex<double> Amplitude;
  typedef bool spin_type;
  typedef State_Ket<spin_type,Amplitude, num_modes,num_bits> State_Ket;
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
  std::ofstream of(params.output_file.c_str());
  if(of){
    params.save(of);
    params.write_header(of);
    for(int i = 0; i <params.N; i++){
      State_Vector v = h.get_psi_lbl();
      std::array<double,2> pop = matrix_recorder.get_spin_pop(v.begin(),v.end());
      params.write_pop_run(of,i+1,i*dt,pop[0],pop[1]);
      h.simple_run(dt);
    }
  }
  of.close();
  
  return 0;
}
