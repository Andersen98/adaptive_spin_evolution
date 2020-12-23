#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "configuration.hpp"
#include "hamiltonian/hamiltonian.hpp"
#include <vector>
#include <string>
#include <complex>
#include <memory>

using namespace std;
namespace py = pybind11;


using state_ket = State_Ket<NUM_MODES,NUM_BITS>;
using state_vector = vector<state_ket>;

PYBIND11_MAKE_OPAQUE(state_ket);
PYBIND11_MAKE_OPAQUE(state_vector);
PYBIND11_MAKE_OPAQUE(hamiltonian);
PYBIND11_MODULE(pyket,m){
  m.doc() = "adaptive spin evolution ported to python";
  m.def("num_modes",[](){return int( NUM_MODES);},"Number of allowed modes in a state");
  m.def("num_bits", [](){return int(NUM_BITS);},"Number of bits that can store a level");
  py::class_<state_vector>(m,"StateVector")
    .def(py::init<>())
    .def("__iter__",[](state_vector &v){
      return py::make_iterator(v.begin(),v.end());
    },py::keep_alive<0,1>())
    .def("__len__",[](const state_vector &v){return v.size();})
    .def("clear", &state_vector::clear)
    .def("push_back", [](state_vector &v, const state_ket &k)mutable{v.push_back(k);} );
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
    .def("get_amp",[](const state_ket &k){return k.amp;})
    .def("set_amp",[](state_ket &k,std::complex<double> a)mutable{k.amp = a;});
  py::class_<hamiltonian>(m,"H")
    .def(py::init([](const std::string &json_arg){
      return std::unique_ptr<hamiltonian>(new hamiltonian(json_arg));
    }))
    .def("run_step",&hamiltonian::run_step)
    .def("get_state_vector", &hamiltonian::get_state_vector)
    .def("get_spin_pop",&hamiltonian::get_spin_pop)
    .def("get_emitter_cavity_prob",&hamiltonian::get_emitter_cavity_prob)
    .def("get_mode_pop",&hamiltonian::get_modeLbl_quanta_pop)
    .def("reset_state",&hamiltonian::reset_state)
    .def("run_evolve",&hamiltonian::run_evolve)
    .def("set_epsilon",&hamiltonian::set_epsilon);
  
}
