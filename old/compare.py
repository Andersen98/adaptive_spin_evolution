import pyket
import qutip
import matplotlib.pyplot as plt
 
import numpy as np
fock_N1 = 40
fock_N2 = 40
#basic enery parameters
w0 = 3
v0 = 4
g0 = .2
w1 = 4
v1 = 5
g1 = .07

#step size
dt = .01



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
dt = .01
N = min(fock_N1,fock_N2)-1
up_list_qu = np.empty(N)
up_list_py = np.empty(N)

t = np.empty(N)
for i in range(min(fock_N1,fock_N2)-1):
    h.run_grow()

for i in range(min(fock_N1,fock_N2)-1):
    print('--------------ITERATION ' +str(i) + '--------------------')    
    h.blas_evolve(complex(0,-dt))
    qinit_state = qinit_state +complex(0,-dt)*H*qinit_state 
    up,down = h.get_blas_spin_pop()
    up_list_py[i] = up
    up_list_qu[i] = abs(qinit_state.unit().ptrace(0)[0][0][0])
    t[i] = dt*i
 
h.reset()
qinit_state = qutip.tensor(spin,vac0,vac1).unit()
for i in range(min(fock_N1,fock_N2)-1):
    h.run_grow()

for i in range(min(fock_N1,fock_N2)-1):
    print('--------------ITERATION ' +str(i) + '--------------------')    
    h.blas_evolve(complex(0,-dt))
    qinit_state = qinit_state +complex(0,-dt)*H*qinit_state 
    up,down = h.get_blas_spin_pop()
    up_list_py[i] = up
    up_list_qu[i] = abs(qinit_state.unit().ptrace(0)[0][0][0])
    t[i] = dt*i

 
h.reset_with_state(pinit_state)
qinit_state = qutip.tensor(spin,vac0,vac1).unit()
for i in range(min(fock_N1,fock_N2)-1):
    h.run_grow()

for i in range(min(fock_N1,fock_N2)-1):
    print('--------------ITERATION ' +str(i) + '--------------------')    
    h.blas_evolve(complex(0,-dt))
    qinit_state = qinit_state +complex(0,-dt)*H*qinit_state 
    up,down = h.get_blas_spin_pop()
    up_list_py[i] = up
    up_list_qu[i] = abs(qinit_state.unit().ptrace(0)[0][0][0])
    t[i] = dt*i
print(h.get_size())

    
fig, ax = plt.subplots(2,1,figsize=(9,4))
ax[0].plot(t,up_list_py,label="pyket_blas_evolve")
ax[0].legend()
ax[1].plot(t,up_list_qu,label="qutip")
ax[1].legend()
plt.show()
