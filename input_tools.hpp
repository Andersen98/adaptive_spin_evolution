#ifndef ADAPTIVE_INPUT_PARAMETERS
#define ADAPTIVE_INPUT_PARAMETERS 1


#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>

struct param_vals{
  std::string atom_path;
  std::string mode_path;
  double energy_cutoff;
  double energy_unit; //specifies a unit of energy defined in eV
  std::string config_file_path;
  int max_occupation;
  //time related values
  double t0;
  double tf;
  double largest_frequency;
  int N;
  //file read in
  std::vector<double> atom_levels;
  std::vector<double> mode_energies;
  std::vector<double> mode_couplings;
};

bool get_params(param_vals &params,int argc, char * argv[]);

struct Spin_Params{

  double energy;
  
};

#endif //ADAPTIVE_INPUT_PARAMETERS 
