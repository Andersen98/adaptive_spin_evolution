import pyket
import qutip

 
import numpy as np
fock_N1 = 2

#basic enery parameters
w0 = 3
v0 = 4
g0 = .2

#step size
dt = .01

#--------qutip code----------
#basic operators

s_eye = qutip.operators.identity(2)
n_eye1 = qutip.operators.identity(fock_N1)


oz = qutip.tensor(qutip.operators.sigmaz(),n_eye1)
oplus = qutip.tensor(qutip.operators.sigmap(),n_eye1)
ox  = qutip.tensor(qutip.operators.sigmax(),n_eye1)
a0 = qutip.tensor(s_eye,qutip.operators.destroy(fock_N1))



#make initial state
vac0 = qutip.states.basis(fock_N1)
spin = qutip.states.basis(2)
qinit_state = qutip.tensor(spin,vac0).unit()



#hamiltonian
H = 0.5*w0*oz + v0*a0.dag()*a0 + g0*ox*(a0.dag()+a0) 



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
    "mode_energies":[v0],
    "mode_couplings":[g0],
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
dt = .1
for i in range(fock_N1-1):
    print('--------------ITERATION ' +str(i) + '--------------------')    
    
    print(h.get_matrix())
    print(H)
    h.run_step(complex(1,0))
    qinit_state = qinit_state +complex(1,0)*H*qinit_state 
    h.run_grow()
    print(h.get_state())
    print(qinit_state)
    print(qinit_state.unit())
