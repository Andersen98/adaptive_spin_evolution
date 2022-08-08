#include "pyket_config.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "configuration.hpp"
#include "hamiltonian.hpp"
#include "input_tools.hpp"
#include <vector>
#include <string>
#include <sstream>
#include <complex>
#include <memory>
#include <algorithm>

PYBIND11_MAKE_OPAQUE(std::vector<int>);
PYBIND11_MAKE_OPAQUE(std::vector<double>);





using namespace std;
namespace py = pybind11;

using state_ket = State_Ket<NUM_MODES,NUM_BITS>;
using state_vector = std::vector<state_ket>;


/* include definitions for module initialization 
defined in main_*****.cpp `
*/
void init_hamiltonian(py::module_ &);
void init_state_ket(py::module_ &);
void init_state_vector(py::module_ &);
void init_params(py::module_ &);


PYBIND11_MODULE(PYKET,m){
  
  /* initialize separately defined module components */
  init_state_ket(m);
  init_state_vector(m);
  init_hamiltonian(m);
  init_params(m);

  m.doc() = "adaptive spin evolution ported to python";
  m.def("num_modes",[](){return int( NUM_MODES);},"Number of allowed modes in a state");
  m.def("num_bits", [](){return int(NUM_BITS);},"Number of bits that can store a level");
  m.def("load_json_str", &load_json_str);
  m.def("sort_state_by_lbl",[](state_vector &v)
	mutable{sort(v.begin(),v.end(),
		     [](auto &it1,const auto &it2){return it1<it2;});
  });
  m.def("sort_state_by_amp",[](state_vector &v)
	mutable{sort(v.begin(),v.end(),
		     [](auto &it1,const auto &it2){return norm(it1.amp)>norm(it2.amp);});
  });
 
 
  py::class_<std::vector<int>>(m, "IntVector")
    .def(py::init<>())
    .def("clear", &std::vector<int>::clear)
    .def("push_back", [](std::vector<int> &v, const int i)mutable{v.push_back(i);} )
    .def("__len__", [](const std::vector<int> &v) { return v.size(); })
    .def("__iter__", [](std::vector<int> &v) {
      return py::make_iterator(v.begin(), v.end());
    }, py::keep_alive<0, 1>());
  
  py::class_<std::vector<double>>(m, "DoubleVector")
    .def(py::init<>())
    .def("clear", &std::vector<double>::clear)
    .def("push_back", [](std::vector<double> &v, const double i)mutable{v.push_back(i);} )
    .def("__len__", [](const std::vector<double> &v) { return v.size(); })
    .def("__iter__", [](std::vector<double> &v) {
      return py::make_iterator(v.begin(), v.end());
    }, py::keep_alive<0, 1>())
    .def("__repr__",[](std::vector<double> &v){
      string result = "[ ";
      for(auto d: v){
	result = result + to_string(d) + ", ";
      }
      result.pop_back();
      result.pop_back();
      result += " ]";
      return(py::str(result));
    });
  

}
