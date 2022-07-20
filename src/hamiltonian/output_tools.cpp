#include "output_tools.hpp"
#include "input_tools.hpp"

#include <boost/format.hpp>

#include <rapidjson/document.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/ostreamwrapper.h>
#include <fstream>



#include <array>
#include <algorithm>
#include <complex>
using namespace std;
using boost::format;
using namespace rapidjson;
namespace adaptive{


  void write_two_vs_time_header(std::ofstream &o, const param_vals &p, string q1, string q2){
    o <<format("PSTART~%|1$-d| DSTART~%|2$-d|\n") %5 %7; 
    format lnBrk("%|88T=|\n");
    o << lnBrk;
    o << format("%|1$=88s|\n")%"RUN BEGIN";
    o << lnBrk;
    o << format("Step(i/%|4$-d|)%|13T~|%|1$-s|%|38T~|%|2$-s|%|63T~|%|3$-s|%|88T~|\n") % "Time" % q1 % q2%p.N;
    o << lnBrk;   
  }

  void write_two_vs_time(std::ofstream &o, const param_vals &p, int id, double time,
				double q1, double q2){
     format
      popFmt("%|1$-d|/%|2$-d|%|13t|%|3$-1.14e|%|38t|%|4$-1.14e|%|63t|%|5$-1.14e|\n");

    o << popFmt % id % p.N % time %q1 %q2; 
    
  }

  
  void write_spin_population_header(std::ofstream &o, const param_vals &p){
  o <<format("PSTART~%|1$-d| DSTART~%|2$-d|\n") %5 %7; 
  format lnBrk("%|88T=|\n");
  o << lnBrk;
  o << format("%|1$=88s|\n")%"RUN BEGIN";
  o << lnBrk;
  o << format("Step(i/%|4$-d|)%|13T~|%|1$-s|%|38T~|%|2$-s|%|63T~|%|3$-s|%|88T~|\n") % "Time" % "Spin Up Pop." % "Spin Down Pop."%p.N;
  o << lnBrk;

  }

  void write_spin_population_run(std::ofstream &o,const param_vals& p, int id, double time,
				 double up, double down){
    format
      popFmt("%|1$-d|/%|2$-d|%|13t|%|3$-1.14e|%|38t|%|4$-1.14e|%|63t|%|5$-1.14e|\n");

    o << popFmt % id % p.N % time %up %down; 
  }

  void write_mode_pop_header(std::ofstream &o, const param_vals &p){
  o <<format("PSTART~%|1$-d| DSTART~%|2$-d|\n") %5 %7; 
  format lnBrk("%|104T=|\n");
  o << lnBrk;
  o << format("%|1$=104s|\n")%"RUN BEGIN";
  o << lnBrk;
  o << format("Step(i/%|5$-d|)%|13T~|%|1$-s|%|38T~|%|2$-s|%|52T~|%|3$-s|%|81T~|%|4$-s|%|98T~|\n") % "Time" % "Mode Label" % "<n>" % "Pop" %p.N;
  o << lnBrk;

  }


  
  void write_mode_pop(std::ofstream &o,const param_vals &p,int id, double time,
				  std::array<std::tuple<int,double,double>,NUM_MODES>
				  mode_lbl_quanta_pop){
    format
      mode("%|6$-d|/%|5$-d|%|13t|%|1$-1.10e|%|38t|%|2$-d|%|52t|%|3$-1.17e|%|81t|%|4$-1.17e|%|104t|\n");
    
    for(int i =0; i < NUM_MODES;i++){
      int lbl = std::get<0>(mode_lbl_quanta_pop[i]);
      double quanta = std::get<1>(mode_lbl_quanta_pop[i]);
      double pop = std::get<2>(mode_lbl_quanta_pop[i]);
      
      o<< mode  %time  % lbl % quanta % pop %p.N %id;
    }  
  }



  void write_stats_header(std::ofstream &o, const param_vals &p){

    format lnBrk("%|88T=|\n");
    o << lnBrk;
    o << format("%|1$=88s|\n")%"STATS BEGIN";
    o << lnBrk;
    o << format("Step(i/%|4$-d|)%|13T~|%|1$-s|%|38T~|%|2$-s|%|63T~|%|3$-s|%|88T~|\n") % "Time" % "Config. Size" % "Exceeded modes"%p.N;
    o << lnBrk; 
  }

  void write_stats(std::ofstream &o, const param_vals &p, int id, double time,
		   int config_space, std::array<int,NUM_MODES> exceeded){
    string exStr = "";
    
    for(int i = 0; i < NUM_MODES; i++){
      int x = exceeded[i];
      if(x>-1){
      	exStr += std::to_string(x)+"/";
      }
    }
    
    format popFmt("%|1$-d|/%|2$-d|%|13t|%|3$-1.14e|%|38t|%|4$-1.14e|%|63t|%|5$-1.14e|\n");
    o << popFmt % id % p.N % time %config_space %exStr; 
  
  }


  void write_state(std::string out_path, const hamiltonian &h){

    const array<int,NUM_MODES> &new2old = h.get_new2old();
    const vector<State_Ket<NUM_MODES,NUM_BITS>> & state_vector_ = h.get_state_vector();
    vector<State_Ket<NUM_MODES,NUM_BITS>> state_vector (state_vector_);
    sort(state_vector.begin(),state_vector.end(),[](auto &it1,auto&it2){return it1.idx < it2.idx;});
    
    //construct document
    Document d;
    d.SetObject();

    auto & allocator = d.GetAllocator();
  

    //construct state_vector array
    Value state_vector_v(kArrayType);
    for(int i = 0; i < int(state_vector.size()); i ++){
      Value re(real(state_vector[i].amp));
      Value im(imag(state_vector[i].amp));
      Value spin(state_vector[i].spin);
      Value active_modes(kArrayType);
      for(int j = 0; j < NUM_MODES;j++){
	int level = state_vector[i].get_mode(j);
	if(level > 0){
	  Value el(kObjectType);
	  el.AddMember("mode",new2old[j], allocator);
	  el.AddMember("level", level,allocator);

	  active_modes.PushBack(el,allocator);
	  
	}
      }

      Value state_ket_v(kObjectType);
      state_ket_v.AddMember("order_found",state_vector[i].idx,allocator);
      state_ket_v.AddMember("re",re,allocator);
      state_ket_v.AddMember("im",im,allocator);
      state_ket_v.AddMember("spin",spin,allocator);
      state_ket_v.AddMember("active_modes",active_modes,allocator);
      
      state_vector_v.PushBack(state_ket_v,allocator);
      
    }

    //add state_vector array to document
    d.AddMember("size",int(state_vector.size()),allocator);
    d.AddMember("state_vector",state_vector_v,allocator);

    //write document to stream
    
    std::ofstream ofs(out_path.c_str());
    OStreamWrapper osw(ofs);

    PrettyWriter<OStreamWrapper> writer(osw);
    d.Accept(writer);

    ofs.close();
    
  }
  
}
