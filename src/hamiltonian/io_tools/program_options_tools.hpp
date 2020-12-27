#ifndef ADAPTIVE_PROGRAM_OPTIONS_TOOLS
#define ADAPTIVE_PROGRAM_OPTIONS_TOOLS

#include <boost/program_options.hpp>
#include "input_tools.hpp"
bool get_params(param_vals &params,int argc, char * argv[]);

struct Spin_Params{

  double energy;
  
};




#endif // ADAPTIVE_PROGRAM_OPTIONS_TOOLS
