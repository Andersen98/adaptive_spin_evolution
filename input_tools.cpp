#include "input_tools.hpp"

namespace pt = boost::property_tree;
namespace po = boost::program_options;
using boost::format;
using namespace std;
using boost::io::group;

void param_vals::save(std::ofstream &o){
  //empty property tree
  pt::ptree tree;

  tree.put("run_info.run_id",run_id);
  tree.put("run_info.output_directory",output_directory);


  //put the simple values into the tree
  tree.put("energy_cutoff",energy_cutoff);
  tree.put("t0",t0);
  tree.put("tf",tf);
  tree.put("energy_unit",energy_unit);
  tree.put("bits_per_mode",NUM_BITS );
  tree.put("dt", dt);

  vector<double>::iterator result;
  result = std::max_element(mode_couplings.begin(),mode_couplings.end());
  int max_idx = std::distance(mode_couplings.begin(),result);

  tree.put("cavity_mode.energy",mode_energies[max_idx]);
  tree.put("cavity_mode.couplings",mode_couplings[max_idx]);
    
  for_each(atom_levels.begin(),atom_levels.end(),
	   [&](double j){tree.add("atom_params.energy_level",j);});

  pt::ptree mode_values;
  for(int i = 0; i < mode_energies.size(); i++){
    pt::ptree child;
    child.put("energy",mode_energies[i]);
    child.put("coupling",mode_couplings[i]);

    mode_values.push_back(std::make_pair("",child));
  }
  pt::ptree mode_info;
  mode_info.put("total_count",mode_energies.size());
  mode_info.put("strongly_coupled", 1);
  mode_info.put("weakly_coupled", mode_energies.size()-1);
  mode_info.put("spectral_density", "random_uniform");
  mode_info.put("spectral_energy", spectral_energy);
  tree.add_child("modes.mode_info",mode_info);
  tree.add_child("modes.mode_values",mode_values);
 
  pt::write_json(o,tree);
}

void param_vals::write_header(ofstream &o){

  format lnBrk("%|88T=|\n");
  o << lnBrk;
  o << format("%|1$=88s|\n")%"RUN BEGIN";
  o << lnBrk;
  o << format("Step(i/%|4$-d|)%|13T~|%|1$-s|%|38T~|%|2$-s|%|63T~|%|3$-s|%|88T~|\n") % "Time" % "Spin Up Pop." % "Spin Down Pop."%N;
  o << lnBrk;

}
//void param_vals::write_pop_run(std::ofstream&, int, double, double, double)
void param_vals::write_pop_run(std::ofstream &o, int id, double time, double up, double down){
  format popFmt("%|1$-d|/%|2$-d|%|13t|%|3$-1.14e|%|38t|%|4$-1.14e|%|63t|%|5$-1.14e|\n");
  o << popFmt % id % N % time %up %down; 
}

void param_vals::write_stats_header(std::ofstream &o){

  format lnBrk("%|88T=|\n");
  o << lnBrk;
  o << format("%|1$=88s|\n")%"STATS BEGIN";
  o << lnBrk;
  o << format("Step(i/%|4$-d|)%|13T~|%|1$-s|%|38T~|%|2$-s|%|63T~|%|3$-s|%|88T~|\n") % "Time" % "Config. Size" % "Exceeded modes"%N;
  o << lnBrk;

}

void param_vals::write_stats(std::ofstream &o, int id, double time,
			     int config_space, std::vector<int> exceeded){
  string exStr = "";
  if(exceeded.size()){
    for(auto &x: exceeded){
      exStr += std::to_string(x)+"/";
    }
  }else{
    exStr = "N/A";
  }
    
  format popFmt("%|1$-d|/%|2$-d|%|13t|%|3$-1.14e|%|38t|%|4$-1.14e|%|63t|%|5$-1.14e|\n");
  o << popFmt % id % N % time %config_space %exStr; 

}


//returns exit status
//if exit status is true, end the program
//if exit status is false, continue with the calculation
bool get_params(param_vals &conf, int argc, char * argv[]){
  
  
  
  po::options_description generic("Generic options");
  generic.add_options()
    ("help,h", "produce help message")
    ("config_file", po::value<string>(&conf.config_file_path), "Path to configuration file. Specify to use configuration file")
    ("print_inputs,P", "print the input values that will be used for calculation")
    ("verbose,v", "Print much more information." )
    ;
  
  po::options_description run_information("Run Information");
  run_information.add_options()
    ("run_id",po::value<int>(&conf.run_id),"ID used to label this run. Specifies name of output file in output directory.")
    ("output_dir",po::value<string>(&conf.output_directory),"Directory to store output file");

  po::options_description input_paths("Input Paths");
  input_paths.add_options()
     ("atom_p,a",po::value<string>(&conf.atom_path)->default_value("atom_energies"),"input file is of format where each atomic energy level is on a new line")
    ("mode_p,m",po::value<string>(&conf.mode_path)->default_value("mode_energies"),"input file to fully specify mode values and coupling terms. Read mode evergy then read coupling strength energy. These are read in pairs. ex \"mode1 g1 mod2 g2\" is a good input");

  po::options_description energy_parameters("Energy Parameters");
  energy_parameters.add_options()
    ("cutoff", po::value<double>(&conf.energy_cutoff)->default_value(0.01),"Energy cutoff for |c_i * H_{ij}|, where c_i is the apmlitude of the ith configuration and H_{ij} is the transition matrix element from the ith configuration to the kth configuration.")
    ("energy_unit, u",po::value<double>(&conf.energy_unit)->default_value(1),"Energy scaling applied to the inputs. Atom and mode energies are divided by this number");
  


  po::options_description time_parameters("Time Parameters");
  time_parameters.add_options()
    ("t0",po::value<double>(&conf.t0)->default_value(0.0),"start time of simulation")
    ("tf", po::value<double>(&conf.tf)->default_value(12.0),"end time of simulation")
    ("dt",po::value<double>(&conf.dt)->default_value(.0001),"Time step taken in solving diffEq. Should be smaller than fastest time scale to avoid aliasing.");

  
    //Add options available to both command line and config file
  po::options_description cmdline_options;
  cmdline_options.add(generic).add(run_information).add(input_paths).add(energy_parameters).add(time_parameters);
  
  po::options_description config_file_options;
  config_file_options.add(run_information).add(input_paths).add(energy_parameters).add(time_parameters);
    
  po::options_description visible("All options");
  visible.add(cmdline_options);
    
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,cmdline_options),vm);
  po::notify(vm);

  if(vm.count("help")){
    cout << visible <<endl;
    return true;
  }
    
  if(vm.count("config_file")){
    bool result = false;
    
    cout << "Searching for config file: " << conf.config_file_path << endl;
    ifstream ifs_param_file(conf.config_file_path.c_str());
    if (!ifs_param_file){
      cout << "Cannot find config file, using command line arguments as input parameters."<<endl;
      result = true;
    }else{
      cout << "Found config file: " << conf.config_file_path << endl;
      cout << "Loading config file: (overwrites command line arguments!)" << endl;
      store(po::parse_config_file(ifs_param_file, config_file_options),vm);
      notify(vm);
    }
    ifs_param_file.close();
    if (result){
      return true;
    }
  }
  
 
  //set time steps
  conf.N = (conf.tf-conf.t0)/conf.dt;
    
  if(vm.count("print_inputs") || vm.count("verbose")){
    cout << "run id: " << conf.run_id << endl;
    cout << "output dir: " << conf.output_directory << endl;
    cout << "output file: " << conf.output_directory + to_string(conf.run_id)+".dat"<< endl;
    cout <<  "atom_file: " << conf.atom_path << endl;
    cout <<  "mode_file: " << conf.mode_path << endl;
    cout <<  "cutoff: " << conf.energy_cutoff << endl;
    cout <<  "energy_unit: " << conf.energy_unit << endl;
    cout << "t0: " << conf.t0 << endl;
    cout << "tf: " << conf.tf << endl;
    cout << "dt: " << conf.dt << endl;
    cout << "N : " << conf.N << endl;
    if(vm.count("print_inputs")){
      return true;
    }
  }

  //read in data for modes and couplings and atomic energies
  cout << "Attempting to load emitter parameters at: " << conf.atom_path <<endl;
  ifstream ifs(conf.atom_path.c_str());
  if(!ifs){
    cerr<< "Cannot open: " << conf.atom_path << endl;
    cerr<< "Aborting run" <<endl;
    ifs.close();
    return false;
  }else{
    cout << "Success\n";
    for(double i; ifs >> i;){
      conf.atom_levels.push_back(i);
      //cout << conf.atom_levels.back();
    }
  }
  ifs.close();


  cout << "Attempting to load mode parameters at: " << conf.mode_path <<endl;
  ifstream ifs_mode(conf.mode_path.c_str());
  if(!ifs_mode){
    cerr<< "Cannot open: " << conf.mode_path << endl;
    cerr<< "Aborting run" <<endl;
    ifs_mode.close();
    return false;
    }else{
    cout << "Success\n";
    for(double i; ifs_mode >> i;){
      conf.mode_energies.push_back(i);
      ifs_mode >> i;
      conf.mode_couplings.push_back(i);
      
      //cout << "mode " << conf.mode_energies.back();
      //cout << "\t coupling " << conf.mode_couplings.back() <<endl;
    }
  }  
  ifs_mode.close();

  
  assert(vm.count("run_id"));
  assert(vm.count("output_dir"));
  assert(NUM_MODES==conf.mode_energies.size());
  assert(NUM_MODES==conf.mode_couplings.size());
  //make sure file system exists
  std::filesystem::create_directory(conf.output_directory);
  
  
  return false; //exit status is false (meaning we should proceed
}
