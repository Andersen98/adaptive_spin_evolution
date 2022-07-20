#ifndef ADAPTIVE_SPIN_CONFIGURATION
#define ADAPTIVE_SPIN_CONFIGURATION 1

#include <complex>
#include <vector>
#include <bitset>
#include <cassert>
#include <iostream>
#include <array>
#include <boost/format.hpp>
#include <string>
#include <utility>



class Configuration
{
public:

  
  
  
  
};

template<const int num_modes,const int num_bits>
class partial_config{
  template<int umode,int ubits>
  friend std::ostream& operator<<( std::ostream&, const partial_config<umode,ubits>& p );
  std::array<long unsigned, (num_modes*num_bits/64) +1> rep; //represented in long ints
  static std::vector<long unsigned> mode_masks;
  static const long unsigned one = 1;
  static const long unsigned zero = 0;
  
public:

  partial_config():rep{}{}
  partial_config(const partial_config& c):rep{0}{
    std::copy(std::begin(c.rep),std::end(c.rep),std::begin(rep));
  }
  


  
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
    std::copy(begin(c.rep),end(c.rep),begin(rep));
    
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
          
};

/*
const int num_modes = 300;
const int num_bits = 4;
const int giant_words = (num_modes*num_bits)/64 + 1;
*/
template<int bits>
std::vector<long unsigned> get_mask(){
  int _num_bits = bits;
  unsigned long mode_tail = 1;
  std::vector<long unsigned> mask(64/_num_bits);
  for(int i = 0; i < (_num_bits-1); i++){
    mode_tail |= mode_tail<<1;
  }
  for(int i = 0; i < 64/_num_bits; i++){
    mask[i] = mode_tail<<(i*_num_bits);
  }
  return mask;
}
template< int modes,int bits>
std::vector<long unsigned> partial_config<modes,bits>::mode_masks = get_mask<bits>();//variable resolution VERY strange

template <int modes,int bits>
std::ostream& operator<<( std::ostream& o, const partial_config<modes, bits>& p ) {
  o << "TODO: impliment this";
  return o;
}







//State_Ket is the main interface. Full config and complex are behind the scenes
template<int num_modes,int num_bits >
class State_Ket{
  template< int num_modes_, int num_bits_>
  friend std::ostream& operator<<( std::ostream&, const State_Ket<num_modes_,num_bits>& p );
public:
  typedef bool Spin_Type ;
  typedef std::complex<double> Amplitude;
  typedef partial_config<num_modes,num_bits> Label;



  constexpr static int null_idx = -1;
  constexpr static int empty_idx = -2;
  
  Spin_Type spin;
  Amplitude amp;
  Label lbl;
  //idx is the absolute labeling
  int idx;
  State_Ket():spin(false),amp(0,0),lbl(),idx(empty_idx){}
  //copy constructor
  State_Ket(const State_Ket&c ):spin(c.spin),amp(c.amp),lbl(c.lbl),idx(c.idx){}
  
  
  
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

  std::string str()const{
    std::string result = "";
    boost::format fmtr("Idx:%|4$-d|  Amp:(%|1$-f|,%|2$-f|)  Spin:%|3$-b| ");
    fmtr % std::real(amp)%std::imag(amp)% spin%idx;
    result = fmtr.str();
    result += lbl.str();
    return(result);
  }
  
};

template<int nm, int nb>
std::ostream& operator<<( std::ostream& o, const State_Ket<nm,nb>& p ){
    
  o << p.str();
  return o;

}





#endif //ADAPTIVE_SPIN_CONFIGURATION
