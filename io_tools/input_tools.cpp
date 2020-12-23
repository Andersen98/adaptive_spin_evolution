#include "input_tools.hpp"
#include <sstream>
#include <istream>
namespace pt = boost::property_tree;
namespace po = boost::program_options;
namespace rj = rapidjson;
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
  tree.put("bits_per_mode",NUM_BITS );
  tree.put("dt", dt);

  vector<double>::iterator result;
  result = std::max_element(mode_couplings.begin(),mode_couplings.end());
  int max_idx = std::distance(mode_couplings.begin(),result);

  tree.put("cavity_mode.energy",mode_energies[max_idx]);
  tree.put("cavity_mode.couplings",mode_couplings[max_idx]);
    
  
  pt::ptree mode_values;
  for(uint i = 0; i < mode_energies.size(); i++){
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
  mode_info.put("energy_spectral_density", energy_spectral_density);
  tree.add_child("modes.mode_info",mode_info);
  tree.add_child("modes.mode_values",mode_values);
 
  pt::write_json(o,tree);
}





bool param_vals::load_json(std::istream &ifs){
  initial_state.clear();

  bool result = false;
  rj::IStreamWrapper isw(ifs);

  rj::Document d;
  d.ParseStream(isw);
  //run info
  const rj::Value &ri = d["run_info"];
  run_id = ri["run_id"].GetInt();
  assert(NUM_MODES==ri["num_modes"].GetInt());
  output_directory = ri["system_paths"]["code_output_dir"].GetString();
  //time info
  const rj::Value &t = d["time_params"];
  t0 = 0;
  tf = t["tf"].GetDouble();
  dt = t["dt"].GetDouble();
  N = tf/dt;
  //initial_state
  const rj::Value &psi = d["initial_state"];
  assert(psi.IsArray());
  for(rj::SizeType i = 0; i < psi.Size(); i++){
    complex<double> amp(psi[i]["re"].GetDouble(),psi[i]["im"].GetDouble());
    bool spin = psi[i]["spin"].GetBool();
    int mode_idx = psi[i]["idx"].GetInt();
    int n = psi[i]["n"].GetInt();

    simple_ket k;
    k.amp = amp;
    k.spin = spin;
    k.mode = mode_idx;
    k.n = n;
    initial_state.push_back(k);
  }

  //----------energy---------
  //#energy cutoff
  const rj::Value &en = d["energy_info"];
  energy_cutoff = en["params"]["cutoff"].GetDouble();
  energy_spectral_density = en["params"]["energy_spectral_density"].GetDouble();
  //#emitter energy
  up_energy = en["energies"]["emitter"]["up"].GetDouble();
  down_energy = en["energies"]["emitter"]["down"].GetDouble();
  //#mode energies
  const rj::Value &m = en["energies"]["modes"];
  assert(m.IsArray());
  assert(m.Size() == NUM_MODES);
  mode_energies.clear();
  mode_energies.resize(NUM_MODES);
  mode_couplings.clear();
  mode_couplings.resize(NUM_MODES);
  for(rj::SizeType i = 0; i < m.Size(); i ++){
    mode_energies[i] = m[i]["w"].GetDouble();
    mode_couplings[i] = m[i]["g"].GetDouble();
  }  
  
  return(result);

}


param_vals load_json_str(const string &s){
  param_vals p;
  istringstream ifs;
  ifs.str(s);
  p.load_json(ifs);
  return p;
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
    ("json_file,j",po::value<string>(&conf.json_path),
     "Path to json file. Load parameters from json file, ignoring parameters from command line")
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
    ("cutoff", po::value<double>(&conf.energy_cutoff)->default_value(0.01),"Energy cutoff for |c_i * H_{ij}|, where c_i is the apmlitude of the ith configuration and H_{ij} is the transition matrix element from the ith configuration to the kth configuration.");


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
    if(vm.count("verbose")){
      cout << "Searching for config file: " << conf.config_file_path << endl;
    }
    ifstream ifs_param_file(conf.config_file_path.c_str());
    if (!ifs_param_file){
      cout << "Cannot find config file, Exiting."<<endl;
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

  if(vm.count("json_file")){
    if(vm.count("verbose")){
      cout << "Attempting to load json configuration at: " << conf.json_path <<endl;
    }
    ifstream ifs(conf.json_path.c_str());
    if(!ifs){
      cerr << "Cannot open json file at " << conf.json_path <<endl;
      cerr << "Aborting run" << endl;
      ifs.close();
      return true;
    }
    else{
      bool failure = conf.load_json(ifs);
      ifs.close();
      return failure;
    }

  }
 
  //read in data for modes and couplings and atomic energies
  vector<double> atom_energies(0);
  if(vm.count("verbose")){
    cout << "Attempting to load emitter parameters at: " << conf.atom_path <<endl;
  }
  ifstream ifs(conf.atom_path.c_str());
  if(!ifs){
    cerr<< "Cannot open: " << conf.atom_path << endl;
    cerr<< "Aborting run" <<endl;
    ifs.close();
    return false;
  }else{
    if(vm.count("verbose")){
      cout << "Success\n";
    }
    for(double i; ifs >> i;){
      atom_energies.push_back(i);

    }
  }
  ifs.close();

  //set up and down energy
  conf.up_energy = atom_energies[0];
  conf.down_energy = atom_energies[1];


  if(vm.count("verbose")){
    cout << "Attempting to load mode parameters at: " << conf.mode_path <<endl;
  }
  
  ifstream ifs_mode(conf.mode_path.c_str());
  if(!ifs_mode){
    cerr<< "Cannot open: " << conf.mode_path << endl;
    cerr<< "Aborting run" <<endl;
    ifs_mode.close();
    return false;
    }else{
    if(vm.count("verbose")){
      cout << "Success\n";
    }
    for(double i; ifs_mode >> i;){
      conf.mode_energies.push_back(i);
      ifs_mode >> i;
      conf.mode_couplings.push_back(i);
      
      //cout << "mode " << conf.mode_energies.back();
      //cout << "\t coupling " << conf.mode_couplings.back() <<endl;
    }
  }  
  ifs_mode.close();

  
  //set time steps
  conf.N = (conf.tf-conf.t0)/conf.dt;
    
  if(vm.count("print_inputs") || vm.count("verbose")){
    cout << "run id: " << conf.run_id << endl;
    cout << "output dir: " << conf.output_directory << endl;
    cout << "output file: " << conf.output_directory + to_string(conf.run_id)+".dat"<< endl;
    cout <<  "atom_file: " << conf.atom_path << endl;
    cout <<  "mode_file: " << conf.mode_path << endl;
    cout << "spin up energy: " << conf.up_energy <<endl;
    cout << "spin down energy: " << conf.down_energy <<endl;
    cout <<  "cutoff: " << conf.energy_cutoff << endl;
    cout << "t0: " << conf.t0 << endl;
    cout << "tf: " << conf.tf << endl;
    cout << "dt: " << conf.dt << endl;
    cout << "N : " << conf.N << endl;
    if(vm.count("print_inputs")){
      return true;
    }
  }

  
  assert(vm.count("run_id"));
  assert(vm.count("output_dir"));
  assert(NUM_MODES==conf.mode_energies.size());
  assert(NUM_MODES==conf.mode_couplings.size());
  //make sure file system exists
  std::filesystem::create_directory(conf.output_directory);
  
  
  return false; //exit status is false (meaning we should proceed
}
