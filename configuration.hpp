#ifndef ADAPTIVE_SPIN_CONFIGURATION
#define ADAPTIVE_SPIN_CONFIGURATION 1

#include <complex>
#include <vector>
#include <iostream>
#include <bitset>
#include <cassert>
#include <iostream>
#include <array>
#include <boost/format.hpp>
#include <string>
#include <utility>
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

//



using namespace std;
template<const int num_modes,const int num_bits>
class partial_config{
  template<int umode,int ubits>
  friend std::ostream& operator<<( std::ostream&, const partial_config<umode,ubits>& p );
  //zero indexing for modes ex; 0,1,2,3...
  //friend std::ostream  &operator<<(std::ostream &os, const partial_config &c); 

  std::array<long unsigned, (num_modes*num_bits/64) +1> rep; //represented in long ints
  static vector<long unsigned> mode_masks;
  static const long unsigned one = 1;
  static const long unsigned zero = 0;
  
public:

  partial_config():rep{0}{
    rep.fill(0);
  }
  partial_config(const partial_config& c):rep{0}{
    std::copy(std::begin(c.rep),std::end(c.rep),std::begin(rep));
  }
  
  constexpr static int max_level_int = (2<<num_bits -1);
  constexpr static int num_modes_int = num_modes;
  constexpr static int num_bits_int = num_bits;
  constexpr static int giant_count = (num_modes*num_bits/64) +1;
  
  

 
  
  bool operator<(const partial_config &d) const{
    //copied from Dice Determinat.h
    
    for (int i= int(giant_count) -1; i >=0; i--){
      if(rep[i] > d.rep[i] /*int compare is swapped from mode compare */) {
	return true;
      }else if(rep[i] < d.rep[i]){
	return false;
      }
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
    copy(begin(c.rep),end(c.rep),begin(rep));
    
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
  for(int i = 0; i < (_num_bits-1); i++){
    mode_tail |= mode_tail<<1;
  }
  for(int i = 0; i < 64/_num_bits; i++){
    mask[i] = mode_tail<<(i*_num_bits);
  }
  return mask;
}
template< int modes,int bits>
vector<long unsigned> partial_config<modes,bits>::mode_masks = get_mask<bits>();//variable resolution VERY strange

template <int modes,int bits>
std::ostream& operator<<( std::ostream& o, const partial_config<modes, bits>& p ) {

  int modes_per_line = 10;
  std::array<std::pair<int,int>,modes> mode_level_array;
  int num_active_modes=0;
  for(int i = 0; i < modes; i++){
    int level = p.get_mode(i);
    if(level > 0){
      num_active_modes++;
      mode_level_array[num_active_modes-1] = make_pair<int,int>(std::forward<int&&>(i),std::forward<int&&>(level)) ;
    }
  }
  int rem =  num_active_modes%modes_per_line;
  int div = num_active_modes/modes_per_line;
  std::string format_str = "%|1$-d|:|2$-d|%|8t|%|3$-d|:|4$-d|%|16t|%|5$-d|:|6$-d|%|24t|%|7$-d|:|8$-d|%|32t|%|9$-d|:|10$-d|%|40t|%|11$-d|:|12$-d|%|48t|%|13$-d|:|14$-d|%|56t|%|15$-d|:|16$-d|%|64t|%|17$-d|:|18$-d|%|72t|%|19$-d|:|20$-d|%|80t|\n";

  boost::format fmtr(format_str);
  boost::format fmtr_single("%|1$-d|:|2$-d|");
  int j = 0;
  
  if(div){
    for(int i = 0; i < modes_per_line*div; i++){
      fmtr % mode_level_array[j].first;
      fmtr % mode_level_array[j].second;
      j++;
    }
    o << fmtr;
  }else if(rem){
      
    for(int i = 0; i < rem; i++){
      o << fmtr_single % mode_level_array[j].first % mode_level_array[j].second << "  ";
    }
  }
  
  return o;
  
}







//State_Ket is the main interface. Full config and complex are behind the scenes
template<int num_modes,int num_bits >
class State_Ket{
  template< int num_modes_, int num_bits_>
  friend std::ostream& operator<<( std::ostream&, const State_Ket<num_modes_,num_bits>& p );
public:

  State_Ket():spin(false),amp(0.0),lbl(),idx(empty_idx){}
  //copy constructor
  State_Ket(const State_Ket&c ):idx(c.idx),spin(c.spin),amp(c.amp),lbl(c.lbl){}
  
  typedef std::complex<double> Amplitude;
  typedef partial_config<num_modes,num_bits> Label;
  typedef bool Spin_Type;

  Spin_Type spin;
  Amplitude amp;
  Label lbl;
  //idx is the absolute labeling
  int idx;
  constexpr static int  num_modes_int = num_modes;
  constexpr static int max_level_int = 1<<num_bits -1;
  constexpr static int null_idx = -1;
  constexpr static int empty_idx = -2;
  void add(State_Ket &c){
    amp += c.amp;
  }
  void subtract(State_Ket &c){
    amp -= c.amp;
  }
  void negate(){
    amp = -amp;
  }
  
  //multiply by scalar
  void multiply (Amplitude s){
    amp *= s;
  }
  
  //conpare label + spin
  bool operator <(const State_Ket &c) const{
    if( spin < c.spin){
      return true;
      
    }else if(spin > c.spin){
      return false;
      
    }else if(lbl < c.lbl){
      return true;
    }else{
      //spin == c.spin, lbl > c.lbl
      //spin == c.spin, lbl == c.lbl
      return false;
    }
  }
    
  

  bool operator==(const State_Ket &c)const{
    return( (spin== c.spin)&& (lbl == c.lbl));
  }

  void operator=(const State_Ket &c){
    amp = c.amp;
    lbl = c.lbl;
    spin = c.spin;
    idx = c.idx;
  }
  void set_mode(int mode,long unsigned level){
    lbl.set_mode(mode,level);
  }
  int get_mode(int mode) const{
    return lbl.get_mode(mode);
  }
  void increment_mode(int mode){
    lbl.increment_mode( mode);
  }
  void decrement_mode(int mode){
    lbl.decrement_mode(mode);
  }

  
};

template<int nm, int nb>
std::ostream& operator<<( std::ostream& o, const State_Ket<nm,nb>& p ){
  boost::format fmtr("Amp:(%|1$-e|,%|2$-e|)  Spin:%|2$-b|\n");
  o << fmtr;
  o << p.lbl;
  return o;

}





#endif //ADAPTIVE_SPIN_CONFIGURATION
