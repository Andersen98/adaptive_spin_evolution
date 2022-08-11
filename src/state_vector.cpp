#include "pyket_config.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <vector>
#include <string>
#include"configuration.hpp"
 
PYBIND11_MAKE_OPAQUE(State_Ket<NUM_MODES,NUM_BITS>);
PYBIND11_MAKE_OPAQUE(std::vector<State_Ket<NUM_MODES,NUM_BITS>>);

namespace py = pybind11;

using state_ket = State_Ket<NUM_MODES,NUM_BITS>;
using state_vector = std::vector<state_ket>;

void init_state_vector(py::module_ &m){
 py::class_<state_vector>(m,"StateVector")
    .def(py::init<>())
    .def("__iter__",[](state_vector &v){
      return py::make_iterator(v.begin(),v.end());
    },py::keep_alive<0,1>())
    .def("__len__",[](const state_vector &v){return v.size();})
    .def("clear", &state_vector::clear)
    .def("push_back", [](state_vector &v, const state_ket &k)mutable{v.push_back(k);} )
    .def("__repr__",[](state_vector &v){
      std::string result = "[\n ";
      for(const auto &d: v){
	result = result + "{" + d.str() + "}\n";
      }
      result += "]";
      return(py::str(result));
    })
    .def("__getitem__",[](const state_vector &v,uint i)
    {return v[i];});
}
 