import numpy as np
import argparse
import matplotlib.pyplot as plt
from enum import Enum,unique


def make_tuples_hermetian(tuple_list):
   rows = [x[0] for x in tuple_list]
   cols = [x[1] for x in tuple_list]
   vals = [x[2] for x in tuple_list]
   for r,c,v in tuple_list:
      if r !=c:
         rows.append(c)
         cols.append(r)
         vals.append(np.conjugate(v))
   return(np.array(rows,int),np.array(cols,int),np.array(vals,complex))

class pTests(Enum):
   p100g = 1
   p2q = 2
   p100q = 3
   p1000 = 4
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
         w1 = 7
         v1 = 6
         g1 = .007
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

         
         tf1 = .003
         tf2 = 12
         dt = .0005
         NPTS = 100
            
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
            up1[i],d = h.get_spin_pop()
            up_qutip1[i],d = get_qtip_pop(qinit_state)
            num1[i] = h.get_size()

            t1[i] = dt*i

         up2 =[0 for x in range(N2)]
         down =[0 for x in range(N2)]
         up_qutip2 = [0 for x in range(N2)]
         num2 =[0 for x in range(N2)]
         up_vac_cav_qutip2= [fock_N1*fock_N2*2 for x in range(N2)]
         
         t2 = [0 for x in range(N2)]
         h.set_zero_except_init()
         from scipy.sparse import csc_matrix
         from scipy.sparse.linalg import expm, expm_multiply
         from scipy.linalg import norm
         from scipy.sparse import diags
         v = h.get_tuples()
         n_size = h.get_size()
         row1 = [x[0] for x in v]
         col1 = [x[1] for x in v]
         data1 = [x[2] for x in v]
         for r,c,dat, in v:
            if r != c:
               row1.append(c)
               col1.append(r)
               data1.append(dat)
         
         istate = np.zeros((n_size,1))
         istate[1] = 1
         A=complex(0,-1)*csc_matrix((data1, (row1, col1)), shape=(n_size,n_size))
         timeSci = np.linspace(0,tf2,NPTS,endpoint=True)
         res =expm_multiply(A,istate, start=0, stop=tf2, num=NPTS, endpoint=True)
         spin_idxs = h.get_spin_idxs()

         v_up = np.zeros(n_size)
         for i in range(n_size):
            v_up[i] = spin_idxs[i][0]
         
         upListSci = np.zeros(len(res))
         up_vac_cav = np.zeros(len(res))
         for i in range(len(res)):
            r = diags(v_up,format="csc").dot(res[i])
            val = norm(r)**2
            upListSci[i] = val
            up_vac_cav[i] = norm(res[i][1])**2

         print(upListSci)
         qinit_state = qutip.tensor(spin,vac0,vac1).unit()
         for i in range(N2):
            print(i)
            qinit_state = (qinit_state +complex(0,-dt)*H*qinit_state ).unit()
            up_qutip2[i],d = get_qtip_pop(qinit_state)
            up_vac_cav_qutip2[i] = abs(qinit_state.overlap(qutip.tensor(spin,vac0,vac1).unit()))**2
            t2[i] = dt*i

         up2 = upListSci
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

         
         ax[0][1].plot(timeSci, up2, label="my code population")
         ax[0][1].plot(t2, up_qutip2, label = "qutip")
         ax[0][1].set_xlabel('Time')
         ax[0][1].set_ylabel('Spin Population during Growth')
         ax[0][1].grid(True)
         ax[0][1].legend()
         
         
         ax[1][1].plot(timeSci, up_vac_cav,label='my code')
         ax[1][1].plot(t2, up_vac_cav_qutip2,label='qutip')
         ax[1][1].set_xlabel('Time')
         ax[1][1].set_ylabel('Spin Up Vacuum Cavity Probability')
         ax[1][1].grid(True)
         ax[1][1].legend()
         fig.savefig('/home/ethan/code/run_output/outp2q.png',dpi=300)



         
      if(bt==pTests.p100q and args.bits==8):
         from scipy.sparse import csc_matrix
         from scipy.sparse.linalg import expm, expm_multiply
         from scipy.linalg import norm
         from scipy.sparse import diags
         import pyket100_8 as pyket
         from numpy.random import default_rng
         rng = default_rng()
         w0 = 3
         v0 = 3
         g0 = 2

         m = np.linspace(0,5,100)
         g = np.linspace(0,.01,100)
         g[0] = g0;
         m[0] = v0;
         tf1 = .001
         tf2 = 1000
         dt = .001
         N2 = 10000
         N1 = int(tf1/dt)

         #grow initial space
         psi_init = pyket.StateVector()
         k = pyket.StateKet()
         k.set_mode(0,1)
         k.set_spin(False)
         k.set_amp(complex(0,0))
         psi_init.push_back(k)
         k = pyket.StateKet()
         k.set_spin(True)
         k.set_amp(complex(1,0))
         psi_init.push_back(k)
         d = {
            "initial_state":psi_init,
            "mode_energies":m,
            "mode_couplings":g,
            "up_energy":.5*w0,
            "down_energy":-.5*w0,
            "energy_cutoff":.00001,
         }
         p = pyket.Params()
         p.load_dict(d)
         h = pyket.H(p)

         
         t1 = np.zeros(N1)
         t2 = np.zeros(N2)
         pop1 = np.zeros((N1,2))
         num1 = np.zeros(N1)
         
         for i in range(N1):
            print(i)
            h.grow()
            h.run_step(complex(0,-dt))
            pop1[i][0],pop1[i][1] = h.get_spin_pop()
            num1[i] = h.get_size()
            t1[i] = i*dt
         h.set_zero_except_init()
         rows,cols,data = make_tuples_hermetian(h.get_tuples())
         data = data*complex(0,-1)
         h_size = h.get_size()
         psi_init2 = np.zeros((h.get_size(),1),complex)
         for ket in psi_init:
            idx = h.get_ket_idx(ket)
            assert(idx > -1)
            psi_init2[idx] = ket.get_amp()
            

         pop2 = np.zeros((N2,2))
         up_vacs = np.zeros(N2)
         H_Sparse = csc_matrix((data,(rows,cols)),shape=(h_size,h_size))
         results = expm_multiply(H_Sparse,psi_init2, start=0,stop=tf2,num=N2,endpoint=True)
         spin_idxs = h.get_spin_idxs()
         up_mask = np.array([x[0] for x in spin_idxs])
         down_mask = np.array([x[1] for x in spin_idxs])
         ket = pyket.StateKet()
         ket.set_spin(True)
         up_vac_idx = h.get_ket_idx(ket)
         if(up_vac_idx < 0):
            print("Error!")
         print(h.get_state())
         t2 = np.linspace(0,tf2,N2)
         for i in range(len(results)):
            r_up = diags(up_mask,format="csc").dot(results[i])
            r_down = diags(down_mask,format="csc").dot(results[i])
            pop2[i] = (norm(r_up)**2,norm(r_down)**2)
            up_vacs[i] = abs(results[i][up_vac_idx])**2
         
         

         
         import matplotlib.ticker as ticker
         fig, ax = plt.subplots(2,2 ,figsize=(12,5))
         fig.subplots_adjust(hspace=.5,wspace=.5)

         ax[0][0].plot(t1,pop1[:,0],marker="*",label="Spin Up")
         ax[0][0].plot(t1,pop1[:,1],marker="1",label="Spin Down")
         ax[0][0].set_xlabel('Time')
         ax[0][0].set_ylabel('Spin Population During Growth')
         ax[0][0].grid(True)
         ax[0][0].set_xlim(left=0)
         ax[0][0].xaxis.set_major_locator(ticker.LinearLocator(5))
         ax[0][0].xaxis.set_minor_locator(ticker.LinearLocator(20))
         ax[0][0].legend()
         
         ax[1][0].plot(t1,num1,marker="*")
         ax[1][0].set_xlabel('Time')
         ax[1][0].set_ylabel('State Space Size')
         ax[1][0].set_xlim(left=0)
         ax[1][0].xaxis.set_major_locator(ticker.LinearLocator(5))
         ax[1][0].xaxis.set_minor_locator(ticker.LinearLocator(20))
         ax[1][0].grid(True)
         
         ax[0][1].plot(t2,pop2[:,0],label='Spin Up')
         ax[0][1].plot(t2,pop2[:,1],label='Spin Down')
         ax[0][1].set_xlabel('Time')
         ax[0][1].set_ylabel('Spin Population During Sparse Exp. Mat.')
         ax[0][1].grid(True)
         ax[0][1].legend()
         
         ax[1][1].plot(t2,up_vacs)
         ax[1][1].set_ylabel('Initial Spin Up State Probability')
         ax[1][1].set_xlabel('Time')
         ax[1][1].grid(True)
         plt.show()

      if(bt==pTests.p1000 and args.bits==8):
         from scipy.sparse import csc_matrix
         from scipy.sparse.linalg import expm, expm_multiply
         from scipy.linalg import norm
         from scipy.sparse import diags
         import pyket1000_8 as pyket
         from numpy.random import default_rng
         
         rng = default_rng()
         w0 = 3
         v0 = 3
         g0 = .3

         m = np.linspace(2,4,pyket.num_modes())
         g = np.linspace(0,.01,pyket.num_modes())
         g[0] = g0;
         m[0] = v0;
         tf1 = .11
         tf2 = 70
         dt = .001
         N2 = 1000
         N1 = int(tf1/dt)

         #grow initial space
         psi_init = pyket.StateVector()
         k = pyket.StateKet()
         k.set_mode(0,1)
         k.set_spin(False)
         k.set_amp(complex(0,0))
         psi_init.push_back(k)
         k = pyket.StateKet()
         k.set_spin(True)
         k.set_amp(complex(1,0))
         psi_init.push_back(k)
         d = {
            "initial_state":psi_init,
            "mode_energies":m,
            "mode_couplings":g,
            "up_energy":.5*w0,
            "down_energy":-.5*w0,
            "energy_cutoff":.00001,
         }
         p = pyket.Params()
         p.load_dict(d)
         h = pyket.H(p)

         
         t1 = np.zeros(N1)
         t2 = np.zeros(N2)
         pop1 = np.zeros((N1,2))
         num1 = np.zeros(N1)
         
         for i in range(N1):
            print(i)
            h.grow()
            h.run_step(complex(0,-dt))
            pop1[i][0],pop1[i][1] = h.get_spin_pop()
            num1[i] = h.get_size()
            t1[i] = i*dt
         h.set_zero_except_init()
         rows,cols,data = make_tuples_hermetian(h.get_tuples())
         data = data*complex(0,-1)
         h_size = h.get_size()
         psi_init2 = np.zeros((h.get_size(),1),complex)
         for ket in psi_init:
            idx = h.get_ket_idx(ket)
            assert(idx > -1)
            psi_init2[idx] = ket.get_amp()
            

         pop2 = np.zeros((N2,2))
         up_vacs = np.zeros((N2,2))
         H_Sparse = csc_matrix((data,(rows,cols)),shape=(h_size,h_size))
         results = expm_multiply(H_Sparse,psi_init2, start=0,stop=tf2,num=N2,endpoint=True)
         spin_idxs = h.get_spin_idxs()
         up_mask = np.array([x[0] for x in spin_idxs])
         down_mask = np.array([x[1] for x in spin_idxs])
         ket = pyket.StateKet()
         ket.set_spin(True)
         up_vac_idx = h.get_ket_idx(ket)
         ket.set_spin(False)
         ket.set_mode(0,1)
         down_vac_idx = h.get_ket_idx(ket)
         if(up_vac_idx < 0):
            print("Error!")
         print(h.get_state())
         t2 = np.linspace(0,tf2,N2)
         for i in range(len(results)):
            r_up = diags(up_mask,format="csc").dot(results[i])
            r_down = diags(down_mask,format="csc").dot(results[i])
            pop2[i] = (norm(r_up)**2,norm(r_down)**2)
            up_vacs[i] = abs(results[i][up_vac_idx])**2,abs(results[i][down_vac_idx])**2
         
         

         
         import matplotlib.ticker as ticker
         fig, ax = plt.subplots(2,2 ,figsize=(12,5))
         fig.subplots_adjust(hspace=.5,wspace=.5)

         ax[0][0].plot(t1,pop1[:,0],marker="*",label="Spin Up")
         ax[0][0].plot(t1,pop1[:,1],marker="1",label="Spin Down")
         ax[0][0].set_xlabel('Time')
         ax[0][0].set_ylabel('Spin Population During Growth')
         ax[0][0].grid(True)
         ax[0][0].set_xlim(left=0)
         ax[0][0].xaxis.set_major_locator(ticker.LinearLocator(5))
         ax[0][0].xaxis.set_minor_locator(ticker.LinearLocator(20))
         ax[0][0].legend()
         
         ax[1][0].plot(t1,num1,marker="*")
         ax[1][0].set_xlabel('Time')
         ax[1][0].set_ylabel('State Space Size')
         ax[1][0].set_xlim(left=0)
         ax[1][0].xaxis.set_major_locator(ticker.LinearLocator(5))
         ax[1][0].xaxis.set_minor_locator(ticker.LinearLocator(20))
         ax[1][0].grid(True)
         
         ax[0][1].plot(t2,pop2[:,0],label='Spin Up')
         ax[0][1].plot(t2,pop2[:,1],label='Spin Down')
         ax[0][1].set_xlabel('Time')
         ax[0][1].set_ylabel('Spin Population During Sparse Exp. Mat.')
         ax[0][1].grid(True)
         ax[0][1].legend()
         
         ax[1][1].plot(t2,up_vacs)
         ax[1][1].set_ylabel('Initial Spin Up State Probability')
         ax[1][1].set_xlabel('Time')
         ax[1][1].grid(True)
         plt.show()


         #calculate emission rate
         
