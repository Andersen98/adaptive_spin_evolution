#include "input_tools.hpp"


namespace po = boost::program_options;
using namespace std;


//returns exit status
//if exit status is true, end the program
//if exit status is false, continue with the calculation
bool get_params(param_vals &conf, int argc, char * argv[]){
  
  
    po::options_description generic("Generic options");
    generic.add_options()
      ("help", "produce help message")
       ("config,f", po::value<string>(&conf.config_file_path)->default_value("parameters.conf"), "name of config file")
      ("print_inputs,P", "prind the input values that will be used for calculation")
      ;
    
    po::options_description config("Configuration");
    config.add_options()
      ("atom_p,a",po::value<string>(&conf.atom_path)->default_value("atom_energies"),"input file is of format where each atomic energy level is on a new line")
      ("mode_p,m",po::value<string>(&conf.mode_path)->default_value("mode_parameters"),"input file to fully specify mode values and coupling terms. Read mode evergy then read coupling strength energy. These are read in pairs. ex \"mode1 g1 mod2 g2\" is a good input")
      ("cutoff,c", po::value<double>(&conf.energy_cutoff)->default_value(0.01),"Energy cutoff for |c_i * H_{ij}|, where c_i is the apmlitude of the ith configuration and H_{ij} is the transition matrix element from the ith configuration to the kth configuration.")
      ("energy_unit, u",po::value<double>(&conf.energy_unit)->default_value(1),"Energy scaling applied to the inputs. Atom and mode energies are divided by this number")
      ("max_occupation,o",po::value<int>(&conf.max_occupation)->default_value(16),"Max occupation level for each mode of the bath")
      ("time_start",po::value<double>(&conf.t0)->default_value(0.0),"start time of simulation")
      ("time_end,t", po::value<double>(&conf.tf)->default_value(12.0),"end time of simulation")
      ("n_steps,N",po::value<int>(&conf.N)->default_value(300),"number of steps simulation takes");

    //Add options available to both command line and config file
    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config);

    po::options_description config_file_options;
    config_file_options.add(config);
    
    po::options_description visible("All options");
    visible.add(generic).add(config);
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc,argv,cmdline_options),vm);
    po::notify(vm);

    if(vm.count("help")){
      cout << visible <<endl;
      return true;
    }

    ifstream ifs_param_file(conf.config_file_path.c_str());
    if (!ifs_param_file){
      cout << "cannot open config file: " << conf.config_file_path << endl;
      return true;
    }else{
      store(po::parse_config_file(ifs_param_file, config_file_options),vm);
      notify(vm);
    }

    if(vm.count("print_inputs")){
      cout <<  "cutoff: " << conf.energy_cutoff << endl;
      cout <<  "max_occupation: " << conf.max_occupation << endl;
      cout <<  "atom_file: " << conf.atom_path << endl;
      cout <<  "mode_file: " << conf.mode_path << endl;
      cout <<  "energy_unit: " << conf.energy_unit << endl;
      cout << "t0: " << conf.t0 << endl;
      cout << "tf: " << conf.tf << endl;
      cout << "N : " << conf.N << endl;
      return true;
    }

    //read in data for modes and couplings and atomic energies
    ifstream ifs(conf.atom_path.c_str());
      
    if(!ifs){
      cout<< "Cannot open: " << conf.atom_path << endl;
    }else{
      for(double i; ifs >> i;)
	conf.atom_levels.push_back(i); 
      ifs.close();
    }

    ifstream ifs_mode(conf.mode_path.c_str());
    if(!ifs_mode){
      cout<< "Cannot open: " << conf.mode_path << endl;
    }else{
      for(double i; ifs_mode >> i;){
	conf.mode_energies.push_back(i);
	ifs_mode >> i;
	conf.mode_couplings.push_back(i);
      }
      ifs.close();
    }  

    return false; //exit status is false (meaning we should proceed
}
