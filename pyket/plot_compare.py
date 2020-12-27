import qutip
import pyket

import matplotlib.pyplot as plt
 
import numpy as np
fock_N1 = 14
fock_N2 = 14
#basic enery parameters
w0 = 3
v0 = 4
g0 = .2
w1 = 4
v1 = 5
g1 = .07

#step size
dt = .01

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



p = pyket.Params()
pinit_state = pyket.StateVector()

ket = pyket.StateKet()
ket.set_spin(True)
ket.set_amp(complex(1,0))
pinit_state.push_back(ket)
ket = pyket.StateKet()
ket.set_mode(0,1)
ket.set_amp(complex(0,0))
pinit_state.push_back(ket)
ket = pyket.StateKet()


d = {
    "run_id":1420,
    "output_directory":"~/code/scratch",
    "initial_state":pinit_state,
    "mode_energies":np.linspace(0,10,1000),
    "mode_couplings":np.linspace(0,.0001,1000),
    "up_energy":.5*w0,
    "down_energy":-.5*w0,
    "energy_cutoff":.00005,
    "dt":.01
    
}
N1 =1
N2 = 10000
N = N1+N2
p = pyket.Params()
p.load_dict(d)
h = pyket.H(p)
dt = 0.001
up_list = np.empty(N)
num_states = np.empty(N)
t = np.empty(N)
print(h.get_state())
print(h.get_matrix())
for i in range(N1):
    h.run_grow()
    h.run_step(complex(0,- dt))
    up,down = h.get_spin_pop()
    num_states[i] = h.get_size()
    t[i] = dt*i
    up_list[i] = up
h.set_zero_except_init()
h.switch_evolve()
i_next = N1
for i in range(N1,N):
    h.blas_evolve(complex(0,-dt))
    if(i%10 == 0):
        up,down= h.get_blas_spin_pop()
        num_states[i_next] = h.get_size()
        up_list[i_next] = up
        t[i_next] = i*dt
        i_next = i_next +1
t= np.resize(t,i_next)
up_list = np.resize(up_list,i_next)
num_states = np.resize(num_states,i_next)
fig, ax = plt.subplots(2,1,figsize=(9,4))  # Create a figure containing a single axes.
ax[0].plot(t,up_list)
ax[1].plot(t,num_states)
plt.show()
print(up_list)
