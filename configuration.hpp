#ifndef ADAPTIVE_SPIN_CONFIGURATION
#define ADAPTIVE_SPIN_CONFIGURATION 1

#include <complex>
#include <vector>
#include <iostream>
#include <bitset>
#include <cassert>
#include <iostream>

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



//State_Ket is the main interface. Full config and complex are behind the scenes
template<typename Spin_Type,typename Amplitude, typename value_type, typename Label >
class State_Ket{
  template<typename sp,typename ap, typename vt, typename lb>
  friend std::ostream& operator<<( std::ostream&, const State_Ket<sp,ap,vt,lb>& p );
public:
  
  typedef Amplitude amplitude_type;
  typedef Label lbl_type;
  typedef Spin_Type spin_type;

  Spin_Type spin;
  Amplitude amp;
  Label lbl;
  
  //abs struct contains operator() that does abs https://stackoverflow.com/questions/1174169/function-passed-as-template-argument
  constexpr static int num_modes = Label::num_modes_int;
  constexpr static int max_level = Label::max_level_int;
    
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
    return ((lbl < c.lbl )&&(spin < c.spin));
  }
  
  //JUST CHECKS TO SEE IF LABELS ARE THE SAME. 
  bool operator==(const State_Ket &c)const{
    return( (spin== c.spin)&& (lbl == c.lbl));
  }
  void operator=(const State_Ket &c){
    amp = c.amp;
    lbl = c.lbl;
    spin = c.spin;
  }
  void set_mode(int mode,long unsigned level){
    lbl.set_mode(mode,level);
  }
  int get_mode(int mode){
    return lbl.get_mode(mode);
  }
  void increment_mode(int mode){
    lbl.increment_mode( mode);
  }
  void decrement_mode(int mode){
    lbl.decrement_mode(mode);
  }

  
};


using namespace std;
template<const int num_modes,const int num_bits,int giant_count>
class partial_config{
  template<int umode,int ubits,int ugiant>
  friend std::ostream& operator<<( std::ostream&, const partial_config<umode,ubits,ugiant>& p );
  //zero indexing for modes ex; 0,1,2,3...
  //friend std::ostream  &operator<<(std::ostream &os, const partial_config &c); 

  long unsigned rep[giant_count] = {0}; //represented in long ints
  static vector<long unsigned> mode_masks;
  static const long unsigned one = 1;
  static const long unsigned zero = 0;
  
public:

  constexpr static int max_level_int = (2<<num_bits -1);
  constexpr static int num_modes_int = num_modes;
  constexpr static int num_bits_int = num_bits;
  constexpr static int giant_count_int = giant_count;
  
  

  partial_config(){}
  
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
template< int modes,int bits, int giant>
vector<long unsigned> partial_config<modes,bits,giant>::mode_masks = get_mask<bits>();//variable resolution VERY strange

template <int modes,int bits,int giant>
std::ostream& operator<<( std::ostream& o, const partial_config<modes, bits, giant>& p ) {

  for(int i = 0; i<giant; i++){
    o <<"Modes: " << (i+1)*64-1 << " - " << i*64  <<" ("
      << bitset<64>(p.rep[i]) << " )" << std::endl ;
  }
  return o;
  
}


template<typename sp,typename ap, typename vt,typename lb>
std::ostream& operator<<( std::ostream& o, const State_Ket<sp,ap,vt,lb>& p ){

  o << "Amplitude : " << p.amp << std::endl;
  o << "Spin(T/F): " <<  p.spin << std::endl;
  o << p.lbl;
  return o;

}
#endif //ADAPTIVE_SPIN_CONFIGURATION
