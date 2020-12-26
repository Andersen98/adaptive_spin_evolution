import pyket
fock_N1 = 3
fock_N2 = 3

#basic enery parameters
w0 = 3
v0 = 4
g0 = .2
w1 = 4
v1 = 5
g1 = .07
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



for i in range(min(fock_N1,fock_N2)-1):
    print('--------------ITERATION ' +str(i) + '--------------------')    
    h.run_grow()
    #h.run_step(complex(1,0))
    print(h.get_state())
  
