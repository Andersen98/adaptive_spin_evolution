import random
import numpy as np
from numpy.random import default_rng
rng = default_rng()

omega_0 = 1.87 #units of eV
cavity_0 = 1.91
g_0 = .141

#we use a sharply peaked lorentzian to describe density of states
T = 1/(100*cavity_0) #this is like putting it in a box of length L (sets freq spaceing)
N = 1000
dw = 1/(N*T) #frequency spacing
q = 10 #quality factor (sharpness of lorentzian from scully dos cavity damped)
#center lorentzian to omega_0
# scale the lorentzian so that at omega_0, g_lorentz[argmin(abs(omega_0 - w))] is g_0

def g_lorentz(w,omega_0,g_0,q,N):
    dw,Nw = (w[1]-w[0],len(w))
    idx_near_w0 = np.argmin(np.abs(w-omega_0))
    near_w0 = w[idx_near_w0] #just consider this the new w0
    unscaled_g = (w/(2*q))/((near_w0-w)**2+(w/(2*q))**2)
    norm = g_0/unscaled_g[idx_near_w0]
    g = norm* unscaled_g 
    return(g)

w = np.linspace(0,1/T,N)
g = g_lorentz(w,omega_0,g_0,q,N)

width = .1
offset = .9999
scale = .5
def w_g_dist(N,w_0,width,scale, offset,g_0):
    w = np.zeros(N)
    for i in range(len(w)-1):
        wDelta = (1- offset*np.exp(-width*(np.abs(w[i]-w_0)**2)))
        w[i+1] = w[i]+wDelta

    closest = np.argmin(np.abs(w-w_0))
    norm = g_0/w[closest]
    g = .01 *norm * w
    g[closest] = g[closest]*100
    return w, g
w ,g = w_g_dist(N,omega_0,width,scale, offset,g_0)
print(w)
print(g)

atom_output = "spin_levels"
mode_output = "mode_and_couplings"
with open(atom_output,"w") as f:
    f.write("0"+"\n")
    f.write(str(omega_0)+ "\n")

with open(mode_output,"w") as f:
    for i in range(len(g)):
        f.write( str(w[i]) +"\t")
        f.write( str(g[i]) +"\n")
        
