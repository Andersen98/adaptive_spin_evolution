import random
import numpy as np
from numpy.random import default_rng
import argparse

parser = argparse.ArgumentParser(description='Generate resevior modes, couplings and atom energy parameters')

parser.add_argument('--w0', type=float, help="Splitting of the 2 level emitter")
parser.add_argument('--g0', type=float, help="cavity mode coupling constant")
parser.add_argument('--v0', type=float, help="cavity mode frequency")
parser.add_argument('--gamma', type=float, help="Spectral intensity of bath")
parser.add_argument('-N', type=int,help="Number of oscillators to use")
parser.add_argument('--emitter_out',type=str,default="emitter_out", help="File output for emitter energies")
parser.add_argument('--bosons_out',type=str,default="bosons_out",help="File output for boson energies and couplings")

args = parser.parse_args()


rng = default_rng()

omega_0 = args.w0 
cavity_0 = args.v0
g_0 = args.g0
gamma = args.gamma #spectral intensity of resevior
N = args.N


width = .1
offset = .9999
scale = .5
def w_g_dist(N,w_0,width,scale, offset,g_0,cavity_0,gamma):
    w = np.zeros(N)
    for i in range(len(w)-1):
        wDelta = (1- offset*np.exp(-width*(np.abs(w[i]-w_0)**2)))
        w[i+1] = w[i]+wDelta

    
    norm = 1/np.sum(w**2,axis=0)
    g = gamma *np.sqrt(norm) * w
    w = np.append(w,cavity_0)
    g = np.append(g,g_0)
    return w, g
w ,g = w_g_dist(N,omega_0,width,scale, offset,g_0,cavity_0,gamma)
#print(w)
#print(g)

atom_output = args.emitter_out
mode_output = args.bosons_out
with open(atom_output,"w") as f:
    f.write("0"+"\n")
    f.write(str(omega_0)+ "\n")

with open(mode_output,"w") as f:
    for i in range(len(g)):
        f.write( str(w[i]) +"\t")
        f.write( str(g[i]) +"\n")
        
