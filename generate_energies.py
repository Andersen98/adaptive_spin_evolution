import random

mode_output = "mode_and_couplings"
atom_output = "spin_levels"
N=100
mode_cap = 100
coupling_cap = 30

#if you wanna do a normal dist
coupling_mean = 30
coupling_std = 8

with open(atom_output,"w") as f:
    f.write("0\n")
    f.write("1\n")

with open(mode_output,"w") as f:
    for i in range(N):
        f.write(str(random.random()*mode_cap) +"\t")
        f.write(str(random.random()*coupling_cap) + "\n")
        
