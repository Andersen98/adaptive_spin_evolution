#include "output_tools.hpp"

using namespace std;
using boost::format;

namespace adaptive{

  void write_spin_population_header(ofstream &o, const param_vals &p){
  o <<format("PSTART~%|1$-d| DSTART~%|2$-d|\n") %5 %7; 
  format lnBrk("%|88T=|\n");
  o << lnBrk;
  o << format("%|1$=88s|\n")%"RUN BEGIN";
  o << lnBrk;
  o << format("Step(i/%|4$-d|)%|13T~|%|1$-s|%|38T~|%|2$-s|%|63T~|%|3$-s|%|88T~|\n") % "Time" % "Spin Up Pop." % "Spin Down Pop."%p.N;
  o << lnBrk;

  }

  void write_spin_population_run(ofstream &o,const param_vals& p, int id, double time,
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
			     int config_space, std::vector<int> exceeded){
    string exStr = "";
    if(exceeded.size()){
      for(auto &x: exceeded){
	exStr += std::to_string(x)+"/";
      }
    }else{
      exStr = "N/A";
    }
    
    format popFmt("%|1$-d|/%|2$-d|%|13t|%|3$-1.14e|%|38t|%|4$-1.14e|%|63t|%|5$-1.14e|\n");
    o << popFmt % id % p.N % time %config_space %exStr; 
  
  }

  
}
