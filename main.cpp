#include <iostream>
#include <bitset>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <fstream>
#include <cmath>
#include <complex>
#include <utility>
#include <boost/program_options.hpp>
#include <string>
#include <cassert>
using namespace std;
namespace po = boost::program_options;


inline int bit_count (long x)
{
  x = (x & 0x5555555555555555ULL) + ((x >> 1) & 0x5555555555555555ULL);
  x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
  x = (x & 0x0F0F0F0F0F0F0F0FULL) + ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL);
  return (x * 0x0101010101010101ULL) >> 56;

  //unsigned int u2=u>>32, u1=u;

  //return __builtin_popcount(u2)+__builtin_popcount(u);
  /*
  u1 = u1
    - ((u1 >> 1) & 033333333333)
    - ((u1 >> 2) & 011111111111);


  u2 = u2
    - ((u2 >> 1) & 033333333333)
    - ((u2 >> 2) & 011111111111);

  return (((u1 + (u1 >> 3))
           & 030707070707) % 63) +
    (((u2 + (u2 >> 3))
      & 030707070707) % 63);
  */
}

struct config_vals{
  string atom_path;
  string mode_path;
  double energy_cutoff;
  double energy_unit; //specifies a unit of energy defined in eV
  string config_file_path;
  int max_occupation;
  //time related values
  double t0;
  double tf;
  int N;
};


template<typename spin_type, const int num_modes,const int num_bits,int giant_count>
class partial_config{
  //represents the mode configuration for a spin value in
  // a spin system coupled to resivoir. 
  //zero indexing for modes ex; 0,1,2,3...
  //friend std::ostream  &operator<<(std::ostream &os, const partial_config &c); 
  
  long unsigned rep[giant_count] = {0}; //represented in long ints
  static vector<long unsigned> mode_masks;
  static const long unsigned one = 1;
  static const long unsigned zero = 0;
  spin_type spin;
public:
  
  spin_type get_spin()const{
    return spin;
  }
  void set_spin(spin_type s){
    spin = s;
  }
  
  bool operator<(const partial_config &d) const{
    //copied from Dice Determinat.h
    for (int i=giant_count -1; i >=0; i--){
      if(rep[i] < d.rep[i]) return true;
      else if (rep[i] > d.rep[i]) return false;
    }
    return false;
  }

  bool operator==(const partial_config &d) const{
    for(int i=giant_count-1; i>=0; i--){
      if(rep[i] != d.rep[i]) return false;
    }
    return true;
  }

  void operator=(const partial_config &c) {
    copy(begin(c),end(c),begin(rep));
    spin = c.spin;
  }

  

  //sets mode value to the input level. Leaves all
  //other modes unchanged.
  void set_mode(int mode,long unsigned level){
    assert(mode < num_modes);
    assert(level < (1<<num_bits));
    int idx = (num_bits*mode)/64;
    int offset = (num_bits*mode)%64;
    int mode_offset = offset/num_bits;
    // "~~~" are arbitrary bits
    //  ->~~~~~|currentmodeval|~~~~~~~~
    // &->11111|00000000000000|11111111
    // =->~~~~~|00000000000000|~~~~~~~~
    // |->00000|destinationval|00000000
    // =->~~~~~|destinationval|~~~~~~~~ //desired result
    rep[idx]= (rep[idx]&(~mode_masks[mode_offset])) | (level<<offset); 
  }

  int get_mode(int mode) const {
    assert(mode < num_modes);
    int idx = (mode*num_bits)/64;
    int offset = (mode*num_bits)%64;
    //  "~~~" arbitrary bits
    //  ~~~~~~|requestedval|~~~~~~
    //  (above)>>offset shift to right by offset
    //  =~~~~~~~~~~~~~|requestedval
    //  &0000000000000|111111111111
    //  =0000000000000|requestedval
    return( (rep[idx]>>offset) & mode_masks[0]);

  }

  //increment mode by 1, leave other modes unchanged
  void increment_mode(int mode){
  
    assert(mode <= num_modes);
    int idx = (num_bits*mode)/64;
    int offset = (num_bits*mode)%64;
    // "~~~" are arbitrary bits
    //  ->~~~~~|currentmodeval|~~~~~~~~
    // +->00000|00000000000001|~~~~~~~~
    // =->~~~~~|incrementedval|~~~~~~~~ //desired result
    rep[idx] += one<<offset;
    
  }

  //increment mode by 1, leave other modes unchanged
  void decrement_mode(int mode){

    assert(mode <= num_modes);
    int idx = (num_bits*mode)/64;
    int offset = (num_bits*mode)%64;
    // "~~~" are arbitrary bits
    //  ->~~~~~|currentmodeval|~~~~~~~~
    //(-)->00000|00000000000001|~~~~~~~~
    // =->~~~~~|incrementedval|~~~~~~~~ //desired result
    rep[idx] -= one<<offset;
    
  }
  

  void print_n_words(int n_giant_words=giant_count) const{
    for(int i = 0; i <n_giant_words;  i++){
      std::cout <<"block: " << i << " " << bitset<64>(rep[i]) << endl;
    }
  }

  void print_members()const{
    std::cout << "num_bits: " << num_bits<<endl;
    std::cout << "num_modes: " << num_modes << endl;
    std::cout << "mask size: " << mode_masks.size() << endl;
    for(int i = 0; i < mode_masks.size(); i++){
      cout << "mode " << i << ": " << std::bitset<64>(mode_masks[i]) << endl;
    }
  }
      
    

  
};

/*
const int num_modes = 300;
const int num_bits = 4;
const int giant_words = (num_modes*num_bits)/64 + 1;
*/
template<int bits>
vector<long unsigned> get_mask(){
  int _num_bits = bits;
  unsigned long mode_tail = 1;
  vector<long unsigned> mask(64/_num_bits);
  for(int i = 0; i < _num_bits; i++){
    mode_tail |= mode_tail<<1;
  }
  for(int i = 0; i < 64/_num_bits; i++){
    mask[i] = mode_tail<<(i*_num_bits);
  }
  return mask;
}
template<typename spin, int modes,int bits, int giant>
vector<long unsigned> partial_config<spin,modes,bits,giant>::mode_masks = get_mask<bits>();//variable resolution VERY strange





int main(int argc, char * argv[]){
  
  //declare config vals, which stores config
  config_vals conf;
  
  try{
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
      return 0;
    }

    ifstream ifs(conf.config_file_path.c_str());
    if (!ifs){
      cout << "cannot open config file: " << conf.config_file_path << endl;
      return 0;
    }else{
      store(po::parse_config_file(ifs, config_file_options),vm);
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
      return 0;
    }
    
  }catch(std::exception& e){
    cout << e.what() << endl;
    return 1;
  }

  vector<double> atom_levels;
  vector<double> mode_energies;
  vector<double> mode_couplings;
  
  //read in modes and energies
  try{
    ifstream ifs(conf.atom_path.c_str());
      
    if(!ifs){
      cout<< "Cannot open: " << conf.atom_path << endl;
    }else{
      for(double i; ifs >> i;)
	atom_levels.push_back(i); 
      ifs.close();
    }

    ifstream ifs_mode(conf.mode_path.c_str());
    if(!ifs_mode){
      cout<< "Cannot open: " << conf.mode_path << endl;
    }else{
      for(double i; ifs_mode >> i;){
	mode_energies.push_back(i);
	ifs_mode >> i;
	mode_couplings.push_back(i);
      }
      ifs.close();
    }  
  }catch(std::exception& e){
    cout << e.what() << endl;
    return 1;
  }


  //make transitions
  vector<int> sorted_mode_idx(mode_energies.size());
  iota(sorted_mode_idx.begin(),sorted_mode_idx.end(),0);

  sort(sorted_mode_idx.begin(),sorted_mode_idx.end(),
       [&x=mode_couplings](auto i1,auto i2){ return x[i1] > x[i2];});

  

  //make full transition list
  //I basically implimented a nested for loop with transform and for each (this is more for my learning than for practicality.
  
  vector<vector<double>> transitions(conf.max_occupation);
  transform(transitions.begin(),transitions.end(),transitions.begin(),
	    [&m=mode_couplings, j=1](vector<double> el) mutable {
	      vector<double> v(m);
	      transform(v.begin(),v.end(),v.begin(), [&j](double g){return sqrt(j)*g;});
	      j++;
	      return v;
	    });
  for(auto &v:transitions){
    for(auto i:sorted_mode_idx){
      //cout << v[i] << " " ;
    }
    //cout << endl;
  }
  for(auto i:sorted_mode_idx){
    //cout << mode_couplings[i] << endl;
    
  }

  const int num_spins = 2;
  const int num_modes = 300;
  const int num_bits = 4;
  const int giant_words = (num_bits*num_modes)/64 +1;
  typedef partial_config<bool,num_modes,num_bits,giant_words> full_config;
  typedef vector<pair<full_config,complex<double>>> state_vec;
  state_vec psi(1);
  psi[0].first.set_spin(true);
  psi[0].second = 1;
  
  double t = conf.t0;
  double dt = (conf.tf-conf.t0)/( double(conf.N));
  while(t < conf.tf){

    state_vec delta_phi();
    for(auto &conf_i:v){

      
      
      
    }

    t +=dt;
  }

  

  return 0;
}
