import numpy as np
import argparse
import matplotlib.pyplot as plt
from enum import Enum,unique



class pTests(Enum):
   p100g = 1
   p2q = 2
    
if __name__ == '__main__':
   parser = argparse.ArgumentParser(description='Test different cases for recursion (repeating behaviour)')
   parser.add_argument('tests', metavar='arg', type=lambda val:pTests[val], nargs='+',
                                        help='tests to choose from')
   parser.add_argument('--bits', dest='bits',type=int,
                        help='number of bits to use', default=8)


   args = parser.parse_args()
       

   for bt in args.tests:
      if(bt==pTests.p100g and args.bits==8):
         
         import pyket100_8 as pyket
         w0 = 2
         v0 = 2
         g0 = .2
         tf = .1
         dt = .01
         k1 = pyket.StateKet()
         k2 = pyket.StateKet()
         k1.set_mode(99,1)
         k1.set_spin(False)
         k2.set_spin(True)
         k1.set_amp(complex(.001,0))
         k2.set_amp(complex(1,0))
         init_state = pyket.StateVector()
         init_state.push_back(k1)
         init_state.push_back(k2)
         d = {
            "run_id":1420,
            "output_directory":"~/code/scratch",
            "initial_state":init_state,
            "mode_energies":np.append(np.linspace(0,10,99), v0),
            "mode_couplings": np.append(np.linspace(.0001,.001,99), g0),
            "up_energy":-.5*w0,
            "down_energy":.5*w0,
            "energy_cutoff":.00001,
            "dt":dt
         }
         N = int(tf/dt)
         p = pyket.Params()
         p.load_dict(d)
         h = pyket.H(p)
         up =[0 for x in range(N)]
         num =[0 for x in range(N)]
         t = [0 for x in range(N)]
         for i in range(N):
            h.grow()
            h.run_step(complex(0,-dt))
            up[i],down = h.get_spin_pop()
            num[i] = h.get_size()
            t[i] = dt*i
         fig,(ax1, ax2) = plt.subplots(2, 1)
         # make a little extra space between the subplots
         fig.subplots_adjust(hspace=0.5)
         
         ax1.plot(t, up)
         ax1.set_xlim(0, t[-1])
         ax1.set_xlabel('Time')
         ax1.set_ylabel('Spin Population')
         ax1.grid(True)
         
         
         ax2.plot(t, num)
         ax2.set_xlim(0, t[-1])
         ax2.set_xlabel('Time')
         ax2.set_ylabel('State_Size')
         ax2.grid(True)
         plt.show()
         
      if(bt==pTests.p2q and args.bits==8):

         import pyket2_8 as pyket
         import qutip
         fock_N1 = 40
         fock_N2 = 40
         
         w0 = 3
         v0 = 4
         g0 = .2
         w1 = 4
         v1 = 5
         g1 = .07
         pinit_state = pyket.StateVector()
         
         ket = pyket.StateKet()
         ket.set_spin(True)
         ket.set_amp(complex('1'))
         pinit_state.push_back(ket)
         ket = pyket.StateKet()
         ket.set_mode(0,1)
         ket.set_amp(complex(0))
         pinit_state.push_back(ket)
         ket = pyket.StateKet()
         
         d = {
            "run_id":1420,
            "output_directory":"~/code/scratch",
            "initial_state":pinit_state,
            "mode_energies":[v0,v1],
            "mode_couplings":[g0,g1],
            "up_energy":.5*w0,
            "down_energy":-.5*w0,
            "energy_cutoff":-1,
            "dt":.01  
         }
         
         p = pyket.Params()
         p.load_dict(d)
         h = pyket.H(p)
         
         #--------qutip code----------
         #basic operators
         
         s_eye = qutip.operators.identity(2)
         n_eye1 = qutip.operators.identity(fock_N1)
         n_eye2 = qutip.operators.identity(fock_N2)
         
         oz = qutip.tensor(qutip.operators.sigmaz(),n_eye1,n_eye2)
         oplus = qutip.tensor(qutip.operators.sigmap(),n_eye1,n_eye2)
         ox  = qutip.tensor(qutip.operators.sigmax(),n_eye1,n_eye2)
         a0 = qutip.tensor(s_eye,qutip.operators.destroy(fock_N1),n_eye2)
         a1 = qutip.tensor(s_eye,n_eye1,qutip.operators.destroy(fock_N2))
         
         
         #make initial state
         vac0 = qutip.states.basis(fock_N1)
         vac1 = qutip.states.basis(fock_N2)
         spin = qutip.states.basis(2)
         qinit_state = qutip.tensor(spin,vac0,vac1).unit()
         

         #hamiltonian
         H = 0.5*w0*oz + v0*a0.dag()*a0 + v1*a1.dag()*a1 + g0*ox*(a0.dag()+a0) + g1*ox*(a1.dag()+a1)

         
         tf1 = .05
         tf2 = 10
         dt = .001
         
            
         N1 = int(tf1/dt)
         N2 = int(tf2/dt)
         def get_qtip_pop( qstate):
            return(abs(qstate.unit().ptrace(0)[0][0][0]),abs(qstate.unit().ptrace(0)[0][0][1]))
            
         up1 =[0 for x in range(N1)]
         up_qutip1 = [0 for x in range(N1)]
         num1 =[0 for x in range(N1)]
         num_qutip1 = [fock_N1*fock_N2*2 for x in range(N1)]
         t1 = [0 for x in range(N1)]
         for i in range(N1):
            print(i)
            h.grow()
            h.run_step(complex(0,-dt))
            qinit_state = qinit_state +complex(0,-dt)*H*qinit_state 
            up1[i],down = h.get_spin_pop()
            up_qutip1[i],down = get_qtip_pop(qinit_state)
            num1[i] = h.get_size()

            t1[i] = dt*i

         up2 =[0 for x in range(N2)]
         up_qutip2 = [0 for x in range(N2)]
         num2 =[0 for x in range(N2)]
         num_qutip2 = [fock_N1*fock_N2*2 for x in range(N2)]
         t2 = [0 for x in range(N2)]
         h.store_vector()
         h.store_matrix()
         print("Switching, Size is {0!s}".format(h.get_size()))   
         qinit_state = qutip.tensor(spin,vac0,vac1).unit()
         for i in range(N2):
            print(i)
            qinit_state = qinit_state +complex(0,-dt)*H*qinit_state 
            up2[i],down = h.evolve_state(dt*i)
            up_qutip2[i],down = get_qtip_pop(qinit_state)
            num2[i] = h.get_size()

            t2[i] = dt*i

         fig, ax = plt.subplots(2,2 ,figsize=(12,5) )
         # make a little extra space between the subplots
         fig.subplots_adjust(hspace=0.5)
         
         ax[0][0].plot(t1, up1, label="my code")
         ax[0][0].plot(t1, up_qutip1, label = "qutip")
         ax[0][0].set_xlim(0, t1[-1])
         ax[0][0].set_xlabel('Time')
         ax[0][0].set_ylabel('Spin Population during Growth')
         ax[0][0].grid(True)
         ax[0][0].legend()
         
         
         ax[1][0].plot(t1, num1,label='my code')
         ax[1][0].plot(t1, num_qutip1,label='qutip')
         ax[1][0].set_xlim(0, t1[-1])
         ax[1][0].set_xlabel('Time')
         ax[1][0].set_ylabel('State_Size')
         ax[1][0].grid(True)
         ax[1][0].legend()

         
         ax[0][1].plot(t2, up2, label="my code")
         ax[0][1].plot(t2, up_qutip2, label = "qutip")
         ax[0][1].set_xlim(0, t2[-1])
         ax[0][1].set_xlabel('Time')
         ax[0][1].set_ylabel('Spin Population during Growth')
         ax[0][1].grid(True)
         ax[0][1].legend()
         
         
         ax[1][1].plot(t2, num2,label='my code')
         ax[1][1].plot(t2, num_qutip2,label='qutip')
         ax[1][1].set_xlim(0, t2[-1])
         ax[1][1].set_xlabel('Time')
         ax[1][1].set_ylabel('State_Size')
         ax[1][1].grid(True)
         ax[1][1].legend()
         fig.savefig('/home/ethan/code/run_output/outp2q.png',dpi=300)



         
