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
    "energy_cutoff":0,
    "dt":.01
    
}

p = pyket.Params()
p.load_dict(d)
h = pyket.H(p)
print(pinit_state)
print(qinit_state)
for i in range(min(fock_N1,fock_N2)-1):
    print('--------------ITERATION ' +str(i) + '--------------------')    
    h.run_grow()
    h.run_step(complex(1,0))
    qinit_state = qinit_state +complex(1,0)*H*qinit_state 
    print(h.get_state())
    print(qinit_state.unit())

