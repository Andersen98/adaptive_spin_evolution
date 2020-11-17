#ifndef ADAPTIVE_INPUT_PARAMETERS
#define ADAPTIVE_INPUT_PARAMETERS 1

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>

//property tree

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <exception>
#include <algorithm>
#include <utility>


//write output
#include <boost/format.hpp>
#include <iomanip>



struct param_vals{

  
  //first save, then append output to file

  
  void save(std::ofstream &o);
  void write_header(std::ofstream &o);
  void write_pop_run(std::ofstream &o, int id, double time, double up, double down);
  //  void wrie_populations(ofstream &o);

  
  //run info
  int run_id;
  std::string output_file;

  //input files
  std::string config_file_path;
  std::string atom_path;
  std::string mode_path;

  //physics (energy)
  double energy_cutoff;
  double energy_unit; //specifies a unit of energy defined in eV
  int max_occupation;
  double spectral_energy;
  std::vector<double> atom_levels;
  std::vector<double> mode_energies;
  std::vector<double> mode_couplings;

  //physics (time)
  double t0;
  double tf;
  double largest_frequency;
  int N;
  
};

bool get_params(param_vals &params,int argc, char * argv[]);

struct Spin_Params{

  double energy;
  
};



#endif //ADAPTIVE_INPUT_PARAMETERS 
