import random
import numpy as np
from numpy.random import default_rng


class mode_generator:

    def __init__(self,w0_,g0_,v0_,spectral_energy_,num_modes_,dt_):

        self.w0=w0_
        self.g0=g0_
        self.v0=v0_
        self.spectral_energy=spectral_energy_
        self.num_modes=num_modes_
        self.dt=dt_
        self.rng=default_rng()
        self.num_modes=num_modes_
        self.fastest_time = 10*dt_

    def w_g_rnd_dist(self):
        N = self.num_modes
        gamma = self.spectral_energy
        v0 = self.v0
        g0 = self.g0
        dt = self.dt
        w = np.linspace(self.fastest_time,0,num=(N-1),endpoint=False)
        g = self.rng.random(N-1)
        norm = np.sum(g**2,axis=0)
        g = gamma*g*np.sqrt(1/norm)
        g = np.append(g,g0)
        w = np.append(w,v0)
        return w,g

