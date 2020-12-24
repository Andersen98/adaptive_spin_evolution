#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "configuration.hpp"
#include "hamiltonian/hamiltonian.hpp"
#include "io_tools/input_tools.hpp"
#include <vector>
#include <string>
#include <complex>
#include <memory>
#include <algorithm>
PYBIND11_MAKE_OPAQUE(std::vector<int>);
PYBIND11_MAKE_OPAQUE(std::vector<double>);
PYBIND11_MAKE_OPAQUE(std::vector<simple_ket>);
PYBIND11_MAKE_OPAQUE(State_Ket<NUM_MODES,NUM_BITS>);
PYBIND11_MAKE_OPAQUE(std::vector<State_Ket<NUM_MODES,NUM_BITS>>);
PYBIND11_MAKE_OPAQUE(hamiltonian);


using namespace std;
namespace py = pybind11;

using state_ket = State_Ket<NUM_MODES,NUM_BITS>;
using state_vector = std::vector<state_ket>;


PYBIND11_MODULE(pyket,m){
  m.doc() = "adaptive spin evolution ported to python";
  m.def("num_modes",[](){return int( NUM_MODES);},"Number of allowed modes in a state");
  m.def("num_bits", [](){return int(NUM_BITS);},"Number of bits that can store a level");
  m.def("load_json_str", &load_json_str);
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
    .def(py::init([](param_vals &p){
      return std::unique_ptr<hamiltonian>(new hamiltonian(p));
    }))
    .def("run_grow_evolve",&hamiltonian::run_grow_evolve)
    .def("run_grow",&hamiltonian::run_grow)
    .def("run_grow",&hamiltonian::run_evolve)
    .def("get_state_vector", &hamiltonian::get_state_vector)
    .def("get_spin_pop",&hamiltonian::get_spin_pop)
    .def("get_emitter_cavity_prob",&hamiltonian::get_emitter_cavity_prob)
    .def("get_mode_pop",&hamiltonian::get_modeLbl_quanta_pop)
    .def("set_zero_except_init",&hamiltonian::set_zero_except_init)
    .def("run_evolve",&hamiltonian::run_evolve)
    .def("set_epsilon",&hamiltonian::set_epsilon)
    .def("par_test_one",&hamiltonian::par_test_one)
    .def("par_test_two",&hamiltonian::par_test_two);
  
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
  

  py::class_<simple_ket>(m,"SimpleKet")
    .def(py::init<>())
    .def_readwrite("n",&simple_ket::n)
    .def_readwrite("mode",&simple_ket::mode)
    .def_readwrite("spin",&simple_ket::spin)
    .def_readwrite("amp",&simple_ket::amp)
    .def("__repr__",[](simple_ket &d){
      string result = "";
      result = "spin(" + to_string(d.spin)+")" + ", mode("+to_string(d.mode)+"):lvl("+
	  to_string(d.n) + ") ";
      return py::str(result);
      
    });
    
  py::class_<std::vector<simple_ket>>(m, "SimpleStateVector")
    .def(py::init<>())
    .def("clear", &std::vector<simple_ket>::clear)
    .def("push_back", [](vector<simple_ket> &v, const simple_ket &k)mutable{v.push_back(k);} )
    .def("__len__", [](const std::vector<simple_ket> &v) { return v.size(); })
    .def("__iter__", [](std::vector<simple_ket> &v) {
      return py::make_iterator(v.begin(), v.end());
    }, py::keep_alive<0, 1>())
    .def("__repr__",[](std::vector<simple_ket> &v){
      string result = "[ ";
      for(auto d: v){
	result = result + "{" + to_string(d.spin) + ", mode("+to_string(d.mode)+"):lvl("+
	  to_string(d.n) + ")}, ";
      }
      result.pop_back();
      result.pop_back();
      result += " ]";
      return(py::str(result));
    });
  


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
	  vector<simple_ket> *v = item.second.cast<vector<simple_ket>*>();
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
