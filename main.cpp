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


using namespace std;


struct Spin_Params{

  double energy;
  
};


template<typename State_Vector,int time_steps>
class Population_Callback{

public:
  double time[time_steps] = {};
  double spin_up[time_steps] = {};
  double spin_down[time_steps] = {};
  
  
  void operator()(const State_Vector &psi, double time, int step_count){
    for(auto &c:psi){
      if(c.spin){
	spin_up[step_count] += real(abs(c.amp)*conj(c.amp));
      }else{
	spin_down[step_count] += real(abs(c.amp)*conj(c.amp));
      }
    }
  }

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
  constexpr double (&abs_ptr)(const std::complex<double>&) = abs;  
  typedef complex<double> Amplitude;
  typedef bool spin_type;
  typedef partial_config< num_modes, num_bits, giant_words> Label;
  typedef State_Ket<spin_type,Amplitude, double, Label, abs_ptr > State_Ket;
  typedef vector<State_Ket> State_Vector;
  typedef vector<double>::iterator Iter_m; //mode iterator
  typedef vector<double>::iterator Iter_c; //coupling iterator
  //setup initial conditions
  State_Vector initial_state(1);
  initial_state[1].spin = true;
  initial_state[1].amp = 1;
  initial_state[1].lbl.set_mode(0,3);

  //make callback to calculate the population
  typedef Population_Callback<State_Vector, time_steps> PC;
  PC callback;
  //make hamiltonian
  hamiltonian<Iter_m,Iter_c,Spin_Params, State_Vector, PC>
    h(params.mode_energies.begin(),params.mode_couplings.begin(),spin_params,
      initial_state,params.energy_cutoff);

  
  for(int i = 0; i <time_steps; i++){
    h.do_run( final_time/(time_steps) , callback );
  }

  return 0;
}
