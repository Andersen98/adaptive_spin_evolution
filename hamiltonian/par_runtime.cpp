#include "hamiltonian.hpp"
#include "configuration_char.hpp"
#include <vector>
#include <algorithm>
#include <iostream>
#include <boost/compute/algorithm/copy.hpp>
#include <boost/compute/algorithm/sort.hpp>
#include <boost/compute/container/vector.hpp>
#include <boost/compute/types/struct.hpp>
#include <boost/compute/algorithm/transform.hpp>
#include <boost/compute/functional/math.hpp>
#include <boost/compute/function.hpp>
#include <boost/compute/types/complex.hpp>

namespace compute = boost::compute;


struct basic_config{
  
  bool spin;
  uint8_t rep[63];
  uint idx;
  
  basic_config():spin(false),rep{1},idx(0){}
  bool operator<(const basic_config c)const{
    			 
    if(c.spin >  spin){
      return true;
    }else if( spin == c.spin){
      for(int i = 0; i < 63; i++){
	if(c.rep[i] >  rep[i]){
	  return true;
	}
	return false;
      }
    }

    return false;
  }

};

  BOOST_COMPUTE_ADAPT_STRUCT(basic_config,basic_config, (spin,rep,idx))

void hamiltonian::par_test_one(){



    BOOST_COMPUTE_FUNCTION(bool, cmp_config, (basic_config b, basic_config a),
		       {
			 
			 if(a.spin > b.spin){
			   return true;
			 }else if(b.spin == a.spin){
			   for(int i = 0; i < 63; i++){
			     if(a.rep[i] > b.rep[i]){
			       return true;
			     }
			     return false;
			   }
			 }else{
			   return false;
			 }
		      });



  
  compute::device device = compute::system::default_device();
  compute::context context(device);
  compute::command_queue queue(context,device);

  
  compute::vector<basic_config> state_vec(context);
  for(int i = 0; i < 100000; i++){
    basic_config k;
    k.rep[i%64] = i;
    state_vec.push_back(k,queue);
  }

  compute::sort(state_vec.begin(),state_vec.end(),cmp_config,queue);

 
  
  
}

void hamiltonian::par_test_two(){
 std::vector<basic_config> state_vec2;
  for(int i = 0; i < 1000000; i++){
    basic_config k;
    k.rep[i%64] = i;
    state_vec2.push_back(k);
  }
  std::sort(state_vec2.begin(),state_vec2.end(),
	    [](const basic_config &it1,const basic_config &it2){return it1<it2;});
  


}
