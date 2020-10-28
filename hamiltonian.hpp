#include <iterator>
#include <algorithm>
#include <utility>
#include <vector>
#include <cmath>
#include <iostream>
#include <complex>
//You need to supply your own state ket and mode and coupling iterators
//State Ket can be an std::vector
template<typename Iter_m, typename Iter_c, typename Spin_params_T,
	 typename State_Vector>
class hamiltonian{
  using State_Ket = typename State_Vector::value_type;
  using Label = typename State_Ket::Label;
  using Amplitude = typename State_Ket::amplitude_type;
  using value_type = typename Amplitude::value_type;
  using mode_val_T =  typename std::iterator_traits<Iter_m>::value_type;
  using g_val_T =  typename std::iterator_traits<Iter_c>::value_type;

  constexpr static int max_level = State_Ket::max_level_int;
  constexpr static int num_modes = State_Ket::num_modes_int; 
 
  g_val_T g[State_Ket::max_level_int][State_Ket::num_modes_int];
  mode_val_T m_g[State_Ket::num_modes_int]; //soted by g
  //sorting modes by dec energy
  mode_val_T m_energy[State_Ket::num_modes_int];

  Amplitude one_i = Amplitude(0,1);

  Spin_params_T spin_params;
  State_Vector psi_lbl;
  State_Vector psi_amp;
  State_Vector psi_delta;
  value_type energy_cutoff;
  
public:
  int level_cap = 0;

  hamiltonian():g{{}},m_g{},m_energy{}{}
  
  void setup(Iter_m m_begin,Iter_c c_begin, Spin_params_T spin_params_,
	      State_Vector initial_state, value_type cutoff){
    energy_cutoff = cutoff;
    psi_lbl.clear();
    psi_amp.clear();
    psi_lbl.push_back(initial_state[0]);
    psi_amp.push_back(initial_state[0]);
    sort(psi_lbl.begin(),psi_lbl.end(),[this](auto &it1,auto&it2){return it1 < it2;});
    //important to sort amp in descending order
    sort(psi_amp.begin(), psi_amp.end(),
	 [this](State_Ket &it1,State_Ket &it2){return abs(it1.amp) > abs(it2.amp);});
        
    spin_params_ = spin_params;
    int idx[num_modes] = {};
    mode_val_T m_copy[num_modes] = {};
    g_val_T g_copy[num_modes] = {};
    transform(begin(idx),end(idx),begin(idx),[j=0](auto &el)mutable{return j++;});
    std::copy(m_begin,m_begin+num_modes,begin(m_energy));
    std::copy(m_begin,m_begin+num_modes,begin(m_copy));
    std::copy(c_begin,c_begin+num_modes,begin(g_copy));
    std::sort(begin(m_energy), end(m_energy),[](auto &it1,auto &it2){return it1>it2;});
    std::sort(begin(idx),end(idx), [&g=g_copy](auto &it1,auto &it2){return g[it1] > g[it2];});

    
    
    for(int i = 0; i < num_modes; i++){
      m_g[i] = m_copy[idx[i]];
      for(int n = 0; n < max_level; n++)
	g[n][i] = sqrt(n)*g_copy[idx[i]];
    }

    
  }


  
  static void normalize_state(State_Vector &p){

    double N = 0;
    for(const State_Ket &k: p){
      N += std::real(std::conj(k.amp)*k.amp);
    }
    N = 1/std::sqrt(N);
    
    std::for_each(p.begin(),p.end(),[N](State_Ket & ket){ket.amp*=N;});

  }

  static double state_magnitude( State_Vector &p){
    
    Amplitude mag = 0;
    for(const State_Ket &k: p){

      mag += k.amp*std::conj(k.amp);
    }
    return std::sqrt(std::real(mag));
    
  }

  State_Vector get_psi_lbl(){
    State_Vector v(psi_lbl.size());
    std::copy(psi_lbl.begin(),psi_lbl.end(),v.begin());
    return v;
  }
    

  void simple_run(value_type dt){

    psi_delta.clear();
    bool stop = false;
    for(auto& ket:psi_amp){
      stop = h_dipole(psi_delta,ket,energy_cutoff,dt);
      if(stop)
	break;
    }
    sort_and_merge_by_label(psi_delta,psi_lbl);
    normalize_state(psi_lbl);
    psi_amp.resize(psi_lbl.size());
    copy(psi_lbl.begin(),psi_lbl.end(),psi_amp.begin());
    sort(psi_amp.begin(),psi_amp.end(),
	 [this](auto&it1,auto&it2){return std::abs(it1.amp) > std::abs(it2.amp);});

  }
  
  //require that value_type and value tag match
  template<typename Recorder_Type>
  void do_run(value_type dt, Recorder_Type &p, value_type value_tag){

    psi_delta.clear();
    bool stop = false;
    for(auto& ket: psi_amp){
      stop = h_dipole(psi_delta,ket, energy_cutoff, dt);
      if (stop)
	break;
    }
     sort_and_merge_by_label(psi_delta,psi_lbl);
     normalize_state(psi_lbl);
     psi_amp.resize(psi_lbl.size());
     copy(psi_lbl.begin(),psi_lbl.end(),psi_amp.begin());
     sort(psi_amp.begin(),psi_amp.end(),
	  [this](auto&it1,auto&it2){return std::abs(it1.amp) > std::abs(it2.amp);});
     p(psi_lbl.begin(),psi_lbl.end(),
       typename std::iterator_traits<typename State_Vector::iterator>::iterator_category());
     
  }

 

  //merges delta into psi
  void sort_and_merge_by_label(State_Vector &delta, State_Vector &psi){

    //do the duplicate problem beforehand
    using sIt = typename State_Vector::iterator;

    //merge the duplicates between delta_psi and ps
    for(int i = 0; i<delta.size(); i++ ){
      pair<sIt,sIt> eq= equal_range(psi.begin(),psi.end(),delta[i],
				       [](auto &it1,auto&it2){return it1 < it2 ;}) ;

      //wanna find 34 in a list
      // 0 23 24 56 (end)
      //eq.first points to 24
      //eq.second points to 56
      
      if(std::distance(eq.first,eq.second) > 1){
	//add amplitude
	(eq.first)->add(delta[i]);
	//remove from delta_psi
	delta.erase(delta.begin()+i);
	//adjust indexing from delete
	i--;
      }
    }

    //sort delta_psi by label
    sort(delta.begin(),delta.end(),
	 [](auto &it1,auto&it2){return it1 < it2 ;});

    //reserve space
    psi.reserve(psi.size() + delta.size());

        //save midpoint for later
    sIt midpoint_marker = psi.end();
    //move the delta elements over to psi using move_iter, which
    //moves elements rather than copying them
    psi.insert(psi.end(),
		   make_move_iterator(delta.begin()),
		   make_move_iterator(delta.end()));
    
    inplace_merge(psi.begin(),midpoint_marker,psi.end(),
		  [](auto &it1,auto&it2){return it1 < it2 ;});// [first,middle) and [middle,last)
    

  }

  //this function modifies the delta vector. I.e. adds on to it. 
  bool h_dipole(State_Vector &delta, State_Ket &c, value_type cutoff, value_type dt){

    //reserve space
    delta.reserve(delta.size()+2*num_modes);
    
    Amplitude amp = c.amp;
    
    if(g[max_level][0]*dt*std::abs(amp) <= cutoff){
      return false;
    }
      
    int j=0;
    Amplitude prefactor =  -one_i*dt*amp;
    while(j<num_modes&&(std::abs(prefactor*g[max_level][j]) > cutoff)){
      
      State_Ket r = c;
      State_Ket l = c;
      r.spin = !r.spin;
      l.spin = !l.spin;
      int level = c.get_mode(j);
      Amplitude factor_r;
      Amplitude factor_l;

      switch(level){
      case(0):
	factor_r = prefactor*g[level+1][j];
	if(std::abs(factor_r)> cutoff){
	  r.multiply(factor_r);
	  r.increment_mode(j);
	  delta.push_back(r);
	}
	break;
      case(max_level):
	factor_l = prefactor*g[level][j];
	if(std::abs(factor_l)>cutoff){
	  l.multiply(factor_l);
	  l.decrement_mode(j);
	  delta.push_back(l);
	  level_cap ++;
	}
	break;
      default:
	factor_r = prefactor*g[level+1][j];
	if(std::abs(factor_r) > cutoff){
	  r.multiply(factor_r);
	  r.increment_mode(j);
	  delta.push_back(r);
	}
	factor_l = prefactor*g[level][j];
	if(std::abs(factor_l)>cutoff){
	  l.multiply(factor_l);
	  l.decrement_mode(j);
	  delta.push_back(l);
	}
	break;
      }//end switch
      
      j++;
    }//end j's

   
    return true;
  }//end func

};
