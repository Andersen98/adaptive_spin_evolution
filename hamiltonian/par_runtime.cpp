#include "hamiltonian.hpp"
#include <boost/compute/core.hpp>
#include <iostream>

namespace compute = boost::compute;
void par_test_one(){

  // get the default device
  compute::device device = compute::system::default_device();

  // print the device's name and platform
  std::cout << "hello from " << device.name();
  std::cout << " (platform: " << device.platform().name() << ")" << std::endl;


}

