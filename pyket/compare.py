import pyket as pyket
import numpy as np


p = pyket.Params()
init_state = pyket.SimpleStateVector()

ket = pyket.SimpleKet()
ket.mode = 0
ket.n = 1
ket.spin = True
ket.amp = complex('1+j')
init_state.push_back(ket)
print(ket)
print(init_state)
d = {
    "run_id":1420,
    "output_directory":"~/code/scratch",
    "initial_state":init_state,
    "mode_energies":[2,3],
    "mode_couplings":[4,5],
    "up_energy":1,
    "down_energy":-1,
    "energy_cutoff":.001,
    "dt":.01
    
}

p.load_dict(d)
print(p.initial_state)
h = pyket.H(p)
h.par_test_two()    


