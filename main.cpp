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
    

  const int time_steps = 2000;
  const double final_time = .0002;
  const double energy_cutoff = .0002;
  const double dt = final_time/time_steps;
  const int num_spins = 2;
  const int num_modes = 1;
  const int num_bits = 8;
  const int max_level = (1<<num_bits) -1;

 
  Spin_Params spin_params;
  spin_params.energy = 1;

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
  //make callback to calculate the population
  typedef Spin_Density_Matrix_Evolution<time_steps> PC;
  PC matrix_recorder(dt);
  //make hamiltonian
  hamiltonian<Iter_m,Iter_c,Spin_Params, State_Vector> h;
  
  h.setup(params.mode_energies.begin(),params.mode_couplings.begin(),spin_params,
      initial_state,energy_cutoff);


 
  for(int i = 0; i <time_steps; i++){
    h.simple_run(dt);
    State_Vector v = h.get_psi_lbl();
    matrix_recorder(v.begin(),v.end(),std::iterator_traits<State_Vector::iterator>::iterator_category());
    
    
  }
  cout << matrix_recorder <<endl;
  cout << h.level_cap;
  return 0;
}
