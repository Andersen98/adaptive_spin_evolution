#include "pyket_config.hpp"
#include<pybind11/pybind11.h>


#include "hamiltonian.hpp"
#include "input_tools.hpp"

PYBIND11_MAKE_OPAQUE(State_Ket<NUM_MODES,NUM_BITS>);
PYBIND11_MAKE_OPAQUE(std::vector<State_Ket<NUM_MODES,NUM_BITS>>);
PYBIND11_MAKE_OPAQUE(hamiltonian);

namespace py = pybind11;


using state_ket = State_Ket<NUM_MODES,NUM_BITS>;
using state_vector = std::vector<state_ket>;

void init_hamiltonian(py::module_ &m){

    py::class_<hamiltonian>(m,"H")
    .def(py::init([](const std::string &json_arg){
      return std::unique_ptr<hamiltonian>(new hamiltonian(json_arg));
    }))
    .def(py::init([](param_vals &p){
      return std::unique_ptr<hamiltonian>(new hamiltonian(p));
    }))
    .def("get_ket_idx",&hamiltonian::get_ket_idx,"returns -1 if it cannot get idx.")
    .def("store_matrix",&hamiltonian::store_matrix,py::call_guard<py::gil_scoped_release>())
    .def("store_vector",&hamiltonian::store_vector,py::call_guard<py::gil_scoped_release>())
    .def("evolve_state",&hamiltonian::evolve_state,py::call_guard<py::gil_scoped_release>())
    .def("evolve_step",&hamiltonian::evolve_step,py::call_guard<py::gil_scoped_release>())
    .def("get_tuples",&hamiltonian::get_tuples)
    .def("get_spin_idxs",&hamiltonian::get_spin_idxs)
    .def("reset",&hamiltonian::reset)
    .def("reset_with_state",&hamiltonian::reset_with_state)
    .def("run_step", &hamiltonian::run_step)
    .def("grow",&hamiltonian::run_grow)
    .def("get_state", &hamiltonian::get_state)
    .def("get_size", &hamiltonian::get_psi_size)
    .def("get_spin_pop",&hamiltonian::get_spin_pop)
    .def("get_emitter_cavity_prob",&hamiltonian::get_emitter_cavity_prob)
    .def("get_mode_pop",&hamiltonian::get_modeLbl_quanta_pop)
    .def("set_zero_except_init",&hamiltonian::set_zero_except_init)
    .def("set_epsilon",&hamiltonian::set_epsilon);
}