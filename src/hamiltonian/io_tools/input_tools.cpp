#include "input_tools.hpp"
#include <sstream>
#include <istream>
namespace rj = rapidjson;
using boost::format;
using namespace std;
using boost::io::group;






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
    state_ket k;
    k.amp = amp;
    k.spin = spin;
    k.idx = state_ket::empty_idx;
    assert(psi[i]["modes"].IsArray());
    for(rj::SizeType j = 0; j < psi[i]["modes"].Size(); i++){
      int mode_idx = psi[i]["modes"][j]["idx"].GetInt();
      int n = psi[i]["modes"][j]["n"].GetInt();
      k.set_mode(mode_idx,n);
    }
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
		

