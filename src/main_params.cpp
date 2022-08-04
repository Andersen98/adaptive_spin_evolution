#include "pyket_config.hpp"
#include<pybind11/pybind11.h>
#include<vector>
#include"configuration.hpp"

#include "hamiltonian.hpp"

PYBIND11_MAKE_OPAQUE(State_Ket<NUM_MODES,NUM_BITS>);
PYBIND11_MAKE_OPAQUE(std::vector<State_Ket<NUM_MODES,NUM_BITS>>);
PYBIND11_MAKE_OPAQUE(hamiltonian);

namespace py = pybind11;

using state_ket = State_Ket<NUM_MODES,NUM_BITS>;
using state_vector = std::vector<state_ket>;

void(py::module_ &m){


  py::class_<param_vals>(m, "Params")
    .def(py::init<>())
    .def("load_dict",[](param_vals &p,py::dict dict){
      for(auto item:dict){
	string key = string(py::str(item.first));
	if(key.compare(string{"run_id"})==0){
	  p.run_id = item.second.cast<int>();
	}else if(key.compare(string{"output_directory"})==0){
	  p.output_directory = item.second.cast<string>();
	}else if(key.compare(string{"energy_cutoff"})==0){
	  p.energy_cutoff = item.second.cast<double>();
	}else if(key.compare(string{"energy_spectral_density"})==0){
	  p.energy_spectral_density = item.second.cast<double>();
	}else if(key.compare(string{"mode_couplings"})==0){
	  p.mode_couplings.clear();
	  for(auto val:item.second){
	    p.mode_couplings.push_back(val.cast<double>()); 
	  }
	  assert(p.mode_couplings.size() == NUM_MODES);
	}else if(key.compare(string{"mode_energies"})==0){
	  p.mode_energies.clear();
	
	  for(auto val:item.second){
	    p.mode_energies.push_back(val.cast<double>()); 
	  }
	  assert(p.mode_energies.size() == NUM_MODES);
	}else if(key.compare(string{"up_energy"})==0){
	  p.up_energy = item.second.cast<double>();
	}else if(key.compare(string{"down_energy"})==0){
	  p.down_energy = item.second.cast<double>();
	}else if(key.compare(string{"t0"})==0){
	  p.t0 = item.second.cast<double>();
	}else if(key.compare(string{"tf"})==0){
	  p.tf = item.second.cast<double>();
	}else if(key.compare(string{"dt"})==0){
	  p.dt = item.second.cast<double>();
	}else if(key.compare(string{"N"})==0){
	  p.N = item.second.cast<double>();
	}else if(key.compare(string{"initial_state"})==0){
	  state_vector *v = item.second.cast<state_vector*>();
	  p.initial_state.resize(v->size());
	  copy(v->begin(),v->end(),p.initial_state.begin());
	}else{
	  py::print(py::str("Key not able to be assigned: " + key));
	}

      }

    })
    .def_readwrite("run_id",&param_vals::run_id)
    .def_readwrite("output_directory",&param_vals::output_directory)
    .def_readwrite("energy_cutoff",&param_vals::energy_cutoff)
    .def_readwrite("energy_spectral_density",&param_vals::energy_spectral_density)
    .def_readwrite("mode_couplings",&param_vals::mode_couplings)
    .def_readwrite("mode_energies",&param_vals::mode_energies)
    .def_readwrite("up_energy",&param_vals::up_energy)
    .def_readwrite("down_energy",&param_vals::down_energy)
    .def_readwrite("t0",&param_vals::t0)
    .def_readwrite("tf",&param_vals::tf)
    .def_readwrite("dt",&param_vals::dt)
    .def_readwrite("N",&param_vals::N)
    .def_readwrite("initial_state",&param_vals::initial_state);
    
    
}