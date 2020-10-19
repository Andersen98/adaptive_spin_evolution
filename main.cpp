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


struct Spin_Params{

  double energy;
  
};


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
    params.energy_cutoff;

  //make transitions
  vector<int> s_idx(params.mode_energies.size());
  iota(s_idx.begin(),s_idx.end(),0);

  sort(s_idx.begin(),s_idx.end(),
       [&x=params.mode_couplings](auto i1,auto i2){ return x[i1] > x[i2];});

  //make full transition list 
  vector<vector<double>> g(params.max_occupation);
  transform(g.begin(),g.end(),g.begin(),
	    [&m=params.mode_couplings, j=1](vector<double> el) mutable {
	      vector<double> v(m);
	      transform(v.begin(),v.end(),v.begin(), [&j](double g){return sqrt(j)*g;});
	      j++;
	      return v;
	    });
    

  const int time_steps = 20;
  const double final_time = 2.0;
  const double energy_cutoff = 20;
  const double dt = time_steps/final_time;
  const int num_spins = 2;
  const int num_modes = 100;
  const int num_bits = 4;
  const int max_level = (1<<num_bits) -1;
  const int giant_words = (num_bits*num_modes)/64 +1;
  assert(g[0].size() == num_modes);
  Spin_Params spin_params;
  spin_params.energy = 1;

  //Complex Number Tools
  //also typedefs
  typedef complex<double> Amplitude;
  typedef bool spin_type;
  typedef partial_config< num_modes, num_bits, giant_words> Label;
  typedef State_Ket<spin_type,Amplitude, double, Label> State_Ket;
  typedef vector<State_Ket> State_Vector;
  typedef vector<double>::iterator Iter_m; //mode iterator
  typedef vector<double>::iterator Iter_c; //coupling iterator
  //setup initial conditions
  vector<Label> labelVec = vector<Label>(3);
  labelVec[0].set_mode(2,2);
  cout << labelVec[0].get_mode(0);
  State_Ket initial_config;
  initial_config.spin = true;
  State_Vector initial_state;
  initial_state.push_back(initial_config);
  initial_state[0].spin = true;
  initial_state[0].amp = 1;
  initial_state[0].lbl.set_mode(0,3);
  cout<<"hello"<<std::endl;
  //make callback to calculate the population
  typedef Spin_Density_Matrix_Evolution<time_steps> PC;
  PC matrix_recorder(dt);
  //make hamiltonian
  hamiltonian<Iter_m,Iter_c,Spin_Params, State_Vector>
    h(params.mode_energies.begin(),params.mode_couplings.begin(),spin_params,
      initial_state,energy_cutoff);


 
  for(int i = 0; i <5; i++){
    h.do_run(dt,matrix_recorder,PC::value_tag );
  }
  cout << matrix_recorder;
  return 0;
}
