#define BOOST_TEST_MODULE dataset_example
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <../configuration.hpp>
#include <../hamiltonian.hpp>
//input utilities
#include <../input_tools.hpp>
#include <string>
#include <iterator>
#include <algorithm>
#include <cmath>

namespace data = boost::unit_test::data;
namespace utf = boost::unit_test;


const int modes = 500;
const int bits = 4;
typedef partial_config<modes,bits> p_conf;
BOOST_TEST_DECORATOR(* utf::label("partial_state_comparison"))
BOOST_DATA_TEST_CASE(partial_config1 /*suite name*/,
		     data::xrange(500)*data::xrange(500)^data::random(1,15),
		     m1,m2,level)
{

  p_conf c1;
  p_conf c2;
  c1.set_mode(m1,level);
  c2.set_mode(m2,level);
  BOOST_TEST((c1 < c2)==(m1 > m2));
}

typedef State_Ket<bool, complex<double>, modes,bits> Ket;
BOOST_TEST_DECORATOR(* utf::label("ket_comparison"))
BOOST_DATA_TEST_CASE(State_Ket_1 /*suite name*/,
		     data::xrange(500)*data::xrange(500)^data::random(1,15)^data::random(0,1)^data::random(0,1),
		     m1,m2,level,spin1,spin2)
{

  Ket k1;
  Ket k2;
  k1.spin = spin1;
  k1.set_mode(m1,level);
  
  k2.spin = spin2;
  k2.set_mode(m2,level);
  //(k1 < k2)==(m1>m2)&&(spin1 < spin2));
  bool kcheck = k1<k2;
  
  BOOST_TEST(kcheck == ((m1>m2)&&(spin1 < spin2)));
}




struct F {
public:
  
  const int num_modes=100;
  const int bits = 4;
  typedef  vector<double>::iterator vIt;
  typedef State_Ket<bool,complex<double>,100,4> State_Ket_T;
  typedef vector<State_Ket_T> State_Vector;
  typedef hamiltonian <vIt,vIt,Spin_Params, State_Vector> hamiltonian_T;
  
  
  F(){
    BOOST_TEST_MESSAGE( "setup fixture" );
    //read in params
    int argc = 1;
    char argv[] = {'p','r','\0'};
    char * argptr[] = {argv};
    bool exit =  get_params(params,argc, argptr);
    BOOST_TEST_REQUIRE(!exit);
    
    Spin_Params spin_params;
    spin_params.energy = 1;

    double energy_cutoff = .1;

    State_Vector initial_state;
    State_Ket_T ket;
    ket.amp=1;
    initial_state.push_back(ket);
  }
	
  ~F()
  { BOOST_TEST_MESSAGE( "teardown fixture" ); }
  void setup(){ 
    
  }
  void teardown(){
  
  }

  param_vals params;
  hamiltonian_T h();
  
  
};

BOOST_TEST_DECORATOR(*utf::tolerance(0.00001) *utf::label("operators"))
BOOST_DATA_TEST_CASE_F(F, vector_normalize,
		       data::xrange(300)^data::random()^data::random()^data::random(1,500),
		       idx,re,im,mul)
{

  State_Vector v(300);
  for( int i = 1; i < 300; i++){
    v[i].amp += 1;
    v[i].amp *=mul;
  }
  BOOST_TEST_MESSAGE( v[idx]);
  v[idx].amp =0 ;
  
  F::hamiltonian_T::normalize_state(v);
  double norm = F::hamiltonian_T::state_magnitude(v);
  
  BOOST_TEST( norm == 1 );

}




struct FSMALL {
public:
  
  const int num_modes=2;
  const int bits = 4;
  typedef  vector<double>::iterator vIt;
  typedef State_Ket<bool,complex<double>,2,4> State_Ket_T;
  typedef vector<State_Ket_T> State_Vector_T;
  typedef hamiltonian <vIt,vIt,Spin_Params, State_Vector_T> hamiltonian_T;

  param_vals params;
  State_Vector_T delta;
  State_Ket_T ket;
  hamiltonian_T h;
  
  FSMALL(){
    BOOST_TEST_MESSAGE( "setup fixture" );
    //read in params
    int argc = 1;
    char argv[] = {'p','r','\0'};
    char * argptr[] = {argv};
    bool exit =  get_params(params,argc, argptr);
    BOOST_TEST_REQUIRE(!exit);
    
    Spin_Params spin_params;
    spin_params.energy = 1;

    double energy_cutoff = .00000000001;

    State_Vector_T initial_state(1);
    initial_state[0].amp=1;

    h.setup(params.mode_energies.begin(),params.mode_couplings.begin(),
		    spin_params,initial_state,energy_cutoff);

        
    h.h_dipole(delta,ket,-1,1);
    h.simple_run(1);
    

  }
	
  ~FSMALL()
  { BOOST_TEST_MESSAGE( "teardown fixture" ); }
  void setup(){ 
    
  }
  void teardown(){
  
  }

  
  
};


BOOST_TEST_DECORATOR(*utf::tolerance(0.00001) *utf::label("operators"))
BOOST_DATA_TEST_CASE_F(FSMALL,single_mode,
		       data::xrange(1),
		       r)
{

  State_Vector_T my_v = delta;
  
  BOOST_TEST(my_v[r].get_mode(0)==1);

}


const int mode_one_expected[] = {0,1,0};
const int mode_two_expected[] = {0,0,1};

BOOST_TEST_DECORATOR(*utf::tolerance(0.00001) *utf::label("auxillary"))
BOOST_DATA_TEST_CASE_F(FSMALL, single_run,
		       data::xrange(3)^data::make(mode_one_expected)^data::make(mode_two_expected),
		       idx,m1,m2)
{

  State_Vector_T v = h.get_psi_lbl();
  BOOST_TEST_MESSAGE(v[idx]);
  BOOST_REQUIRE(v.size() == 3);
  BOOST_TEST(v[idx].get_mode(0) == m1);
  BOOST_TEST(v[idx].get_mode(1) == m2);
    
}
		       
