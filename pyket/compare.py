import pyket
import qutip
import matplotlib.pyplot as plt
 
import numpy as np
fock_N1 = 6
fock_N2 = 6
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
    "energy_cutoff":-1,
    "dt":.01
    
}

p = pyket.Params()
p.load_dict(d)
h = pyket.H(p)
print(pinit_state)
print(qinit_state)
dt = .1
for i in range(min(fock_N1,fock_N2)-1):
    h.run_grow()
h.switch_evolve()
for i in range(min(fock_N1,fock_N2)-1):
    print('--------------ITERATION ' +str(i) + '--------------------')    
    h.blas_evolve(complex(0,-dt))
    qinit_state = qinit_state +complex(0,-dt)*H*qinit_state 
    print(h.get_blas_spin_pop())
    print(qinit_state.unit().ptrace(0))

