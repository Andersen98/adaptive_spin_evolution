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
  m.def("sort_state_by_lbl",[](state_vector &v)
	mutable{sort(v.begin(),v.end(),
		     [](auto &it1,const auto &it2){return it1<it2;});
  });
  m.def("sort_state_by_amp",[](state_vector &v)
	mutable{sort(v.begin(),v.end(),
		     [](auto &it1,const auto &it2){return norm(it1.amp)>norm(it2.amp);});
  });
 
  py::class_<state_vector>(m,"StateVector")
    .def(py::init<>())
    .def("__iter__",[](state_vector &v){
      return py::make_iterator(v.begin(),v.end());
    },py::keep_alive<0,1>())
    .def("__len__",[](const state_vector &v){return v.size();})
    .def("clear", &state_vector::clear)
    .def("push_back", [](state_vector &v, const state_ket &k)mutable{v.push_back(k);} )
    .def("__repr__",[](state_vector &v){
      string result = "[\n ";
      for(const auto &d: v){
	result = result + "{" + d.str() + "}\n";
      }
      result += "]";
      return(py::str(result));
    })
    .def("__getitem__",[](const state_vector &v,uint i)
    {return v[i];});
  
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
  py::class_<hamiltonian>(m,"H")
    .def(py::init([](const std::string &json_arg){
      return std::unique_ptr<hamiltonian>(new hamiltonian(json_arg));
    }))
    .def(py::init([](param_vals &p){
      return std::unique_ptr<hamiltonian>(new hamiltonian(p));
    }))
    .def("run_step", &hamiltonian::run_step)
    .def("run_grow",&hamiltonian::run_grow)
    .def("get_state", &hamiltonian::get_state)
    .def("get_size", &hamiltonian::get_psi_size)
    .def("get_spin_pop",&hamiltonian::get_spin_pop)
    .def("get_emitter_cavity_prob",&hamiltonian::get_emitter_cavity_prob)
    .def("get_mode_pop",&hamiltonian::get_modeLbl_quanta_pop)
    .def("set_zero_except_init",&hamiltonian::set_zero_except_init)
    .def("set_epsilon",&hamiltonian::set_epsilon)
    .def("run_ode_evolve",&hamiltonian::odeint_evolve)
    .def("get_matrix",[](hamiltonian &h){
      // Redirect cout.
      streambuf* oldCoutStreamBuf = cout.rdbuf();
      ostringstream strCout;
      cout.rdbuf( strCout.rdbuf() );

      // This goes to the string stream.
      cout << h.get_matrix() ;

      // Restore old cout.
      cout.rdbuf( oldCoutStreamBuf );
      
      
      return strCout.str();
    })
    .def("switch_evolve",&hamiltonian::switch_evolve)
    .def("blas_evolve",&hamiltonian::blas_evolve)
    .def("get_blas_spin_pop",&hamiltonian::get_blas_spin_pop);
  
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
