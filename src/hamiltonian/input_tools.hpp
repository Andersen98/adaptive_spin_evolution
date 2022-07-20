#ifndef ADAPTIVE_INPUT_PARAMETERS
#define ADAPTIVE_INPUT_PARAMETERS 1
#include <cassert>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <istream>
#include <algorithm>
#include <utility>
#include <iomanip>
#include <filesystem>
#include <tuple>
#include <complex>
#include "configuration.hpp"



#include <boost/format.hpp>


#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/ostreamwrapper.h>
#include <rapidjson/writer.h>

  
struct param_vals{
  typedef State_Ket<NUM_MODES,NUM_BITS> state_ket;
  //rapid json
  bool load_json(std::istream &is);
  
  
  //run info
  int run_id;
  std::string output_directory;

  //input files
  std::string config_file_path;
  std::string json_path;
  std::string atom_path;
  std::string mode_path;

  //physics (energy)
  double energy_cutoff;
  double energy_spectral_density;
  std::vector<double> mode_couplings;
  std::vector<double> mode_energies;
  double up_energy;
  double down_energy;
  
  //physics (time)
  double t0;
  double tf;
  double dt;
  int N;

  //initial state
  std::vector<state_ket> initial_state;
  
};
param_vals load_json_str(const std::string &s);






#endif //ADAPTIVE_INPUT_PARAMETERS 
