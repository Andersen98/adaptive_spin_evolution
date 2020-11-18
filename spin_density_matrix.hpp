#ifndef ADAPTIVE_SPIN_DENSITY_MATRIX
#define ADAPTIVE_SPIN_DENSITY_MATRIX 1
#include <utility>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <array>
#include <complex>
#include <vector>
template<typename value_type>
class Spin_Density_Matrix_Evolution{

  template<typename vt>
  friend std::ostream& operator<<( std::ostream& o, const Spin_Density_Matrix_Evolution<vt>& p );
  std::vector<value_type> spin_up; 
  std::vector<value_type> spin_down;
  int time_step = 0;
  int time_steps;
  double dt;
  

public:
  constexpr static double value_tag = 1.1;


  Spin_Density_Matrix_Evolution(double dt_,int time_steps_): dt(dt_),spin_up(time_steps_),spin_down(time_steps_),time_steps(time_steps_){} 

  template <typename State_Iterator>
  void operator ()(State_Iterator begin, State_Iterator end,
		    std::random_access_iterator_tag){
    
    while(begin !=end){
      if(begin->spin){
	spin_up[time_step] +=std::abs(begin->amp )*std::abs(begin->amp );
      }else{
	spin_down[time_step]+= std::abs(begin->amp )*std::abs(begin->amp );
      }

      begin++;
    }
    time_step++;
  }
  template <typename State_Iterator>
  std::pair<double,double> get_spin_pop(State_Iterator begin, State_Iterator end){
    std::pair<double,double> result=std::make_pair<double,double>(0,0);
    while(begin != end){
      if(begin->spin){
	result.first += std::abs(begin->amp)*std::abs(begin->amp);
      }else{
	result.second += std::abs(begin->amp)*std::abs(begin->amp);
      }
      begin++;
    }
    return result;
  }
};



template<typename vt>
std::ostream& operator<<( std::ostream& o, const Spin_Density_Matrix_Evolution<vt>& p ){
  for(int i = 0; i < p.time_steps; i++){
    o << (i+1)*p.dt << " " << p.spin_up[i] << " " <<  p.spin_down[i] << std::endl; 
  }  
  return o;
}

#endif
