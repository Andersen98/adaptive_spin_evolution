#include "hamiltonian.hpp"
#include <boost/numeric/odeint.hpp>
#include <functional>
using namespace std;
typedef boost::numeric::ublas::vector<complex<double>> state_type;

using namespace boost::numeric::odeint;
using namespace boost::numeric::ublas;


typedef runge_kutta_cash_karp54< state_type > error_stepper_type;

void hamiltonian::odeint_evolve(double t1,double dt){
  
  
  
 
  runge_kutta4< state_type > stepper;
  state_type x(psi_amp.size());
  for(auto &k:psi_amp){
    x[k.idx] = k.amp;

  }
  //  state_type dxdt = complex<double>(0,-1)*prod(x,connection_matrix);
  
  integrate_const(stepper,
		     [&](const state_type &y,state_type &dxdt, const double ){
		       dxdt = complex<double>(0,-1)*prod(y,connection_matrix);},x,0.0,t1,dt);

  for(auto &k:psi_amp){
    k.amp = x[k.idx];
    
  }
}

