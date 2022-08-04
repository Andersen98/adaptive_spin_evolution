#include "pyket_config.hpp"
#include<pybind11/pybind11.h>

#include"configuration.hpp"

PYBIND11_MAKE_OPAQUE(State_Ket<NUM_MODES,NUM_BITS>);
PYBIND11_MAKE_OPAQUE(std::vector<State_Ket<NUM_MODES,NUM_BITS>>);


namespace py = pybind11;

using state_ket = State_Ket<NUM_MODES,NUM_BITS>;

namespace py = pybind11;

void init_state_ket(py::module_ &m){
  
 py::class_<state_ket>(m,"StateKet")
    .def(py::init<>())
    .def("get_mode",&state_ket::get_mode)
    .def("set_mode",&state_ket::set_mode)
    .def("increment_mode",&state_ket::increment_mode)
    .def("decrement_mode",&state_ket::decrement_mode)
    .def("get_spin",[](const state_ket &k){return k.spin;})
    .def("set_spin",[](state_ket &k,bool spin)mutable{ k.spin = spin;})
    .def("__repr__",&state_ket::str)
    .def("__eq__",&state_ket::operator==)
    .def("__lt__",&state_ket::operator<)
    .def("set_amp",[](state_ket &k,int i)mutable{k.idx = i;})
    .def("get_idx",[](const state_ket &k){return k.idx;})
    .def("get_amp",[](const state_ket &k){return k.amp;})
    .def("set_amp",[](state_ket &k,std::complex<double> a)mutable{k.amp = a;});
 

  
}