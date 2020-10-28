#ifndef ADAPTIVE_SPIN_DENSITY_MATRIX
#define ADAPTIVE_SPIN_DENSITY_MATRIX 1
#include <iterator>
#include <algorithm>
#include <iostream>
#include <array>
#include <complex>
template<int time_steps>
class Spin_Density_Matrix_Evolution{

  template<int tSteps>
  friend std::ostream& operator<<( std::ostream& o, const Spin_Density_Matrix_Evolution<tSteps>& p );
  std::array<double, time_steps> spin_up = {};
  std::array<double, time_steps> spin_down = {};
  int time_step = 0;
  double dt;
  

public:
  constexpr static double value_tag = 1.1;


  Spin_Density_Matrix_Evolution(double dt_): dt(dt_){} 

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

};


template<int tSteps>
std::ostream& operator<<( std::ostream& o, const Spin_Density_Matrix_Evolution<tSteps>& p ){
  for(int i = 0; i < tSteps; i++){
    o << (i+1)*p.dt << " " << p.spin_up[i] << " " <<  p.spin_down[i] << std::endl; 
  }  
  return o;
}

#endif
