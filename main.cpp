#include <iostream>
#include <bitset>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <cmath>
#include <complex>
#include <utility>
#include <string>
#include <tuple>
#include "configuration.hpp"
#include "io_tools/input_tools.hpp"
#include "io_tools/output_tools.hpp"
#include "hamiltonian/hamiltonian.hpp"


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
    

  
  
  

  
  
  //make hamiltonian
  hamiltonian h(params);
  const double dt = params.dt;


  //setup output stream
  string output_file = params.output_directory+to_string(params.run_id)+".out";
  string stats_file = params.output_directory+to_string(params.run_id)+".stats";
  string mode_file = params.output_directory+to_string(params.run_id)+".modes";
  string final_state_file = params.output_directory+to_string(params.run_id)+".end_state";
  std::ofstream of(output_file.c_str());
  std::ofstream of_stats(stats_file.c_str());
  std::ofstream of_mode(mode_file.c_str());
  vector<int> exceeded(0);
  if(of&& of_stats){
    adaptive::write_spin_population_header(of,params);
    adaptive::write_stats_header(of_stats,params);
    adaptive::write_mode_pop_header(of_mode,params);
    for(int i = 0; i <params.N; i++){
      std::pair<double,double> pop = h.get_spin_pop();
      adaptive::write_spin_population_run(of,params,i+1,i*dt,pop.first,pop.second);
      adaptive::write_stats(of_stats,params, i+1, i*dt,h.get_psi_size(), exceeded);
      //params.write_mode_pop(of_mode,i+1,i*dt,h.get_modeLbl_quanta_pop());
      h.run_step(dt);
    }
  }

  adaptive::write_state(final_state_file,h);
  of.close();
  of_stats.close();
  of_mode.close();
  return 0;
}
