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
  tree.put("run_info.output_file",output_file);


  //put the simple values into the tree
  tree.put("energy_cutoff",energy_cutoff);
  tree.put("t0",t0);
  tree.put("tf",tf);
  tree.put("energy_unit",energy_unit);
  tree.put("max_occupation", max_occupation);
  tree.put("dt", 1/largest_frequency);

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
      ("n_steps,N",po::value<int>(&conf.N)->default_value(300),"number of steps simulation takes")
      ("largest_frequency",po::value<double>(&conf.largest_frequency)->default_value(100),"this should be greater than your highest frequency in the run to avoid aliasing.")
      ("output_file",po::value<string>(&conf.output_file)->default_value("default_code_output"),"File that will be written out to");
    
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
    ifs_param_file.close();

    ofstream ofs_output_file(conf.output_file.c_str());
    if(!ofs_output_file){
      cout << "cannot open output file: " << conf.output_file <<endl;
      return true;
    }
    ofs_output_file.close();
    
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
      return false;
    }else{
      for(double i; ifs >> i;){
	conf.atom_levels.push_back(i);
	//cout << conf.atom_levels.back();
      }
    }
    ifs.close();

    ifstream ifs_mode(conf.mode_path.c_str());
    if(!ifs_mode){
      cout<< "Cannot open: " << conf.mode_path << endl;
      return false;
    }else{
      for(double i; ifs_mode >> i;){
	conf.mode_energies.push_back(i);
	ifs_mode >> i;
	conf.mode_couplings.push_back(i);

	//cout << "mode " << conf.mode_energies.back();
	//cout << "\t coupling " << conf.mode_couplings.back() <<endl;
      }
    }  
    ifs_mode.close();
    
    return false; //exit status is false (meaning we should proceed
}
