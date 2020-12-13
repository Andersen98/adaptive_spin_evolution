#ifndef ADAPTIVE_OUTPUT_PARAMETERS
#define ADAPTIVE_OUTPUT_PARAMETERS 1


#include <cassert>
#include <vector>
#include <array>
#include <utility>
#include <tuple>
#include <string>
#include <iostream>
#include <fstream>

#include <boost/format.hpp>


#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/ostreamwrapper.h>
#include <rapidjson/writer.h>

#include "input_tools.hpp"

namespace adaptive{

  
  void write_spin_population_header(std::ofstream &o, const param_vals&);
  void write_spin_population_run(std::ofstream &o,const param_vals&, int id, double time, double up, double down);
  void write_stats_header(std::ofstream &o, const param_vals&);
  void write_stats(std::ofstream &o, const param_vals&,int id, double time, int config_space, std::vector<int> exceeded);
  void write_mode_pop_header(std::ofstream &o,const param_vals& );
  void write_mode_pop(std::ofstream &o, const param_vals&, int id, double time,std::array<std::tuple<int,double,double>,NUM_MODES> mode_lbl_quanta_pop);




}


#endif //ADAPTIVE_OUTPUT_PARAMETERS 
