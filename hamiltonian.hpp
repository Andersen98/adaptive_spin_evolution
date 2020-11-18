#include <iterator>
#include <algorithm>
#include <utility>
#include <vector>
#include <cmath>
#include <iostream>
#include <complex>
#include <array>
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
  std::array<int, State_Ket::num_modes_int> energy_to_mode_idx;
  
  Amplitude one_i = Amplitude(0,1);

  Spin_params_T spin_params;
  State_Vector psi_lbl;
  State_Vector psi_amp;
  State_Vector psi_delta;
  value_type energy_cutoff;


  
  
public:
  std::vector<int> exceeded ; //flag for if mode was exceeded
  int level_cap = 0;
  
  
  hamiltonian():g{{}},m_g{},m_energy{},energy_to_mode_idx{},exceeded(){}
  
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
    std::transform(begin(idx),end(idx),begin(idx),[j=0](auto &el)mutable{return j++;});
    std::copy(m_begin,m_begin+num_modes,begin(m_energy));
    std::copy(m_begin,m_begin+num_modes,begin(m_copy));
    std::copy(c_begin,c_begin+num_modes,begin(g_copy));
    std::sort(begin(m_energy), end(m_energy),[](auto &it1,auto &it2){return it1>it2;});
    std::sort(begin(idx),end(idx), [&g_copy](auto &it1,auto &it2){return g_copy[it1] > g_copy[it2];});

    
    
    for(int i = 0; i < num_modes; i++){
      m_g[i] = m_copy[idx[i]];
      for(int n = 0; n < max_level; n++){
	g[n][i] = sqrt(n)*g_copy[idx[i]];
      }
    }
    //idx_base = 0,1,2,3,4...n_modes represents default g ordering
    int idx_base[num_modes] = {};
    std::transform(begin(idx_base),end(idx_base),begin(idx_base),[j=0](auto &el)mutable{return j++;});
    std::sort(begin(idx_base),end(idx_base), [&](auto &it1, auto&it2) {return m_g[it1] > m_g[it2] ;});
    std::copy(begin(idx_base),end(idx_base), begin(energy_to_mode_idx));
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

  State_Vector &get_psi_lbl(){
    State_Vector v(psi_lbl.size());
    return psi_lbl;
  }
    

  bool h_field(State_Vector &delta, State_Ket ket, value_type cutoff, value_type dt){

    bool result = false;
    
    value_type prefactor = dt * std::real(std::conj(ket.amp)*ket.amp);

    if( prefactor * max_level*m_energy[num_modes-1]* prefactor < cutoff){
      result = true;
    }
    
    for(int i = 0; (max_level*m_energy[i]* prefactor > cutoff)&&(i < State_Ket::num_modes_int); i++){
      int mode_idx = energy_to_mode_idx[i];
      int level = ket.get_mode(mode_idx);
      if( level > 0){
	value_type energy = m_energy[i]*level;
	State_Ket c = ket;
	c.amp += -dt*one_i*energy;
	delta.push_back(c);
      }
    }


    return result;

  }

  bool h_atom(State_Vector delta, State_Ket ket, value_type cutoff, value_type dt){
    bool result = false;

    
    value_type prefactor = .5*spin_params.energy* dt * std::real(std::conj(ket.amp)*ket.amp);
    
    if (prefactor  < cutoff){
      result = true;
    }else if(ket.spin){

      State_Ket c = ket;
      c.amp *= - one_i* dt * .5*spin_params.energy;
      delta.push_back(c);
    }else if(!ket.spin){
      State_Ket c = ket;
      c.amp *= one_i * dt * .5 * spin_params.energy;
      delta.push_back(c);
    }
      
    return result;
  }

  void simple_run(value_type dt){
    exceeded.clear();
    psi_delta.clear();
    psi_delta.reserve(2*psi_amp.size() );
        
    bool stop = false;
    for(auto& ket:psi_amp){
      stop = h_dipole(psi_delta,ket,energy_cutoff,dt);
      if(stop)
	break;
    }
    /* stop = false;
    for(auto& ket:psi_amp){
      stop = h_field(psi_delta, ket, energy_cutoff,dt);
      if (stop)
	break;
    }
    stop = false;
    for(auto& ket:psi_amp){
      stop = h_atom(psi_delta, ket, energy_cutoff, dt);
	if (stop)
	  break;
    }*/

    
    sort_and_merge_by_label(psi_delta,psi_lbl);
    normalize_state(psi_lbl);
    psi_amp.resize(psi_lbl.size());
    std::copy(psi_lbl.begin(),psi_lbl.end(),psi_amp.begin());
    std::sort(psi_amp.begin(),psi_amp.end(),
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
    /*std::cout <<"PSI BEGIN-----------------"<<std::endl;
    for(auto z:psi){
      std::cout <<z <<endl;
    }
    std::cout <<"delta BEGIN-----------------"<<std::endl;
    for(auto z:delta){
      std::cout <<z <<endl;
      }*/
    
    
    //do the duplicate problem beforehand
    using sIt = typename State_Vector::iterator;

    int orig = delta.size();
    //merge the duplicates between delta_psi and ps
    for(int i = 0; i<delta.size(); i++ ){
      pair<sIt,sIt> eq= equal_range(psi.begin(),psi.end(),delta[i]) ;

      //wanna find 34 in a list
      // 0 23 24 56 (end)
      //eq.first points to 24
      //eq.second points to 56
      /*
      std::cout <<"delta REF-----------------"<<std::endl;
      std::cout << delta[i] << std::endl;
      while(eq.first < eq.second){
	std::cout <<"EQL RANGE:" <<std::endl;
	std::cout << (*eq.first) <<std::endl;
	eq.first++;
	}*/
	
      if(eq.second == psi.begin()){
	//no elements greater than or equal to
      }else if( eq.second == (++psi.begin())){
	//there are no elements greater than, so we have a match!\
	//(case where entire list is a match)
	(eq.first)->amp += delta[i].amp;
	delta.erase(delta.begin()+i);
	i--;

      }else if((eq.second!=psi.begin()) ){
	//This should be the case where we also have a match

	(eq.first)->amp +=delta[i].amp;
	delta.erase(delta.begin()+i);

	//adjust indexing from delete
	i--;
      }
    }


        
    //sort delta_psi by label
    sort(delta.begin(),delta.end());

    //reserve space
    psi.reserve(psi.size() + delta.size());

    //save midpoint for later
    sIt midpoint_marker = psi.end();
    //move the delta elements over to psi using move_iter, which
    //moves elements rather than copying them
    psi.insert(psi.end(),
		   make_move_iterator(delta.begin()),
		   make_move_iterator(delta.end()));
    
    inplace_merge(psi.begin(),midpoint_marker,psi.end());// [first,middle) and [middle,last)
    
    

  }

  //this function modifies the delta vector. I.e. adds on to it. 
  bool h_dipole(State_Vector &delta, State_Ket &c, value_type cutoff, value_type dt){
    bool stop = false;
    //reserve space
    delta.reserve(delta.size()+2*num_modes);
    
    Amplitude amp = c.amp;
    
    if(g[max_level][0]*dt*std::abs(amp) <= cutoff){
      return false;
    }
      
    int j=0;
    Amplitude prefactor =  -one_i*dt*amp;
    while(j < num_modes &&(std::abs(prefactor*g[max_level][j]) > cutoff)){
      
      State_Ket r = c;
      State_Ket l = c;
      State_Ket d = c;
      int level = c.get_mode(j);
      Amplitude factor_r=0;
      Amplitude factor_l=0;

      switch(level){
      case(0):
	factor_r = prefactor*g[level+1][j];
	if(std::abs(factor_r)> cutoff){
	  r.amp = factor_r;
	  r.spin=!r.spin;
	  r.increment_mode(j);
	  delta.push_back(r);
	  d.amp =  prefactor*(m_g[j]*level +  spin_params.energy*((d.spin)? .5 : -.5));
	  delta.push_back(d);
	}
	break;
      case(max_level):
	
	exceeded.push_back(j);
	factor_l = prefactor*g[level][j];
	if(std::abs(factor_l)>cutoff){
	  l.amp = factor_l * prefactor*(m_g[j]*level + spin_params.energy*((l.spin)? .5:-.5));
	  d.amp =  prefactor*(m_g[j]*level +  spin_params.energy*((r.spin)? .5 : -.5));
	  l.spin = !l.spin;
	  l.decrement_mode(j);
	  delta.push_back(l);
	  d.amp =  prefactor*(m_g[j]*level +  spin_params.energy*((d.spin)? .5 : -.5));
	  delta.push_back(d);
	  
	}
	break;
      default:
	bool passed = false;
	factor_r = prefactor*g[level+1][j];
	if(std::abs(factor_r) > cutoff){
	  r.amp = factor_r + prefactor*(m_g[j]*level+ spin_params.energy*((r.spin)? .5:-.5));
	  r.spin = !r.spin;
	  r.increment_mode(j);
	  delta.push_back(r);
	  passed = true;
	}
	factor_l = prefactor*g[level][j];
	if(std::abs(factor_l)>cutoff){
	  l.amp = factor_l * prefactor*(m_g[j]*level + spin_params.energy* ((r.spin)?.5:-.5));
	  l.spin = !l.spin;
	  l.decrement_mode(j);
	  delta.push_back(l);
	  passed = true;
	}
	if(passed){
	  d.amp =  prefactor*(m_g[j]*level +  spin_params.energy*((d.spin)? .5 : -.5));
	  delta.push_back(d);
	}

	break;
      }//end switch
      
      j++;
    }//end j's

   
    return stop;
  }//end func


};
