#include <iterator>
#include <algorithm>
#include <utility>
#include <vector>
#include <cmath>


//You need to supply your own state ket and mode and coupling iterators
//State Ket can be an std::vector
template<typename Iter_m, typename Iter_c, typename Spin_params_T,
	 typename State_Vector ,typename Callback_Class>
class hamiltonian{
  using State_Ket = typename State_Vector::value_type;
  using Label = typename State_Ket::lbl_type;
  using Amplitude = typename State_Ket::amplitude_type;
  using value_type = typename Amplitude::value_type;
  using mode_val_T =  typename std::iterator_traits<Iter_m>::value_type;
  using g_val_T =  typename std::iterator_traits<Iter_c>::value_type;

  constexpr static int max_level = State_Ket::max_level;
  constexpr static int num_modes = State_Ket::num_modes; 
 
  g_val_T g[State_Ket::max_level][State_Ket::num_modes];
  mode_val_T m_g[State_Ket::num_modes]; //soted by g
  //sorting modes by dec energy
  mode_val_T m_energy[State_Ket::num_modes];

  Amplitude one_i = Amplitude(0,1);
  struct abs_struct{
    value_type operator()(Amplitude z){
      return State_Ket::abs_ptr(z);
    }    
  };

  abs_struct abs;
  Spin_params_T spin_params;
  State_Vector psi_lbl;
  State_Vector psi_amp;
  State_Vector psi_delta;
  value_type energy_cutoff;
  value_type total_time;
  int step_count;

public:
  hamiltonian(Iter_m m_begin,Iter_c c_begin, Spin_params_T spin_params_,
	      State_Vector initial_state, value_type cutoff){
    energy_cutoff = cutoff;
    total_time = 0;
    step_count = 0;
    psi_lbl = State_Vector(initial_state);
    psi_amp = State_Vector(initial_state);
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
      for(int n = 0; n <= max_level; n++)
	g[n][i] = sqrt(n)*g_copy[idx[i]];
    }
  }

  
  void do_run(value_type dt, Callback_Class &callback){

    bool stop = false;
    for(auto& ket: psi_amp){
      stop = h_dipole(psi_delta,ket, energy_cutoff, dt);
      if (stop)
	break;
    }
     sort_and_merge_by_label(psi_delta,psi_lbl);
     copy(psi_lbl.begin(),psi_lbl.end(),psi_amp.begin());
     sort(psi_amp.begin(),psi_amp.end(),
	  [this](auto&it1,auto&it2){return abs(it1.amp) > abs(it2.amp);});
     total_time += dt;
     step_count += 1;
     callback(psi_lbl, total_time, step_count);
      
  }

 

  //merges delta into psi
  void sort_and_merge_by_label(State_Vector &delta, State_Vector &psi){

    //do the duplicate problem beforehand
    using sIt = typename State_Vector::iterator;
    vector<pair<sIt,sIt>> eql_itors;
    //merge the duplicates between delta_psi and ps
    for(int i = 0; i<delta.size(); i++ ){
      pair<sIt,sIt> eql= equal_range(psi.begin(),psi.end(),delta[i],
				       [](auto &it1,auto&it2){return it1 < it2 ;}) ;
      if(eql.first !=eql.second){
	//found a dup
	eql_itors.push_back(eql);
	//add amplitude
	(eql.first)->add(delta[i]);
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

    
    Amplitude amp = c.amp;
    if(g[max_level][0]*dt*abs(amp) <= cutoff){
      return false;
    }
      
    int j=0;
    int level = c.get_mode(j);
    Amplitude prefactor =  -one_i*dt*amp;
    while(abs(prefactor*g[max_level][j]) > cutoff){
      State_Ket r = c;
      State_Ket l = c;
      r.spin = !r.spin;
      l.spin = !l.spin;
      Amplitude factor = prefactor*g[level][j];
      
      if (abs(factor) < cutoff)
	level = -1;
      
      switch(level){
      case(-1):
	//don't add anything
	break;
      case(0):
	r.multiply(factor);
	r.increment_mode(j);
	delta.push_back(r);
	break;
      case(max_level):
	l.multiply(factor);
	l.decrement_mode(j);
	delta.push_back(l);
      default:
	r.multiply(factor);
	r.increment_mode(j);
	l.multiply(factor);
	l.decrement_mode(j);
	delta.push_back(r);
	delta.push_back(l);
	break;
      }//end switch
      
      j++;
      level = c.get_mode(j);
    }//end j's

    return true;
  }//end func

};
