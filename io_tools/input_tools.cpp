#include "input_tools.hpp"
#include <sstream>
#include <istream>
namespace pt = boost::property_tree;
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
		

