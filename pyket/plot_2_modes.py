import matplotlib.pyplot as plt
import numpy as np

import pyket2_8 as pyket

# basic enery parameters
w0 = 1.0
v0 = 1.0
v1 = 0.0
g0 = 0.5
g1 = 0.0

# setup run
dt = 0.01
p = pyket.Params()

ket1 = pyket.StateKet()
ket1.set_spin(True)
ket1.set_amp(complex(1, 0))
ket1.set_mode(0, 0)
ket2 = pyket.StateKet()
ket2.set_spin(True)
ket2.set_mode(0, 1)
ket2.set_amp(complex(0, 0))
ket3 = pyket.StateKet()
ket3.set_spin(False)
ket3.set_amp(complex(0,0))
ket3.set_mode(0,1)
ket4 = pyket.StateKet()
ket4.set_spin(False)
ket4.set_mode(0,0)
ket4.set_amp(complex(0,0))

pinit_state = pyket.StateVector()
pinit_state.push_back(ket1)
pinit_state.push_back(ket2)
pinit_state.push_back(ket3)
pinit_state.push_back(ket4)

d = {
    "run_id": 1420,
    "output_directory": "~/code/scratch",
    "initial_state": pinit_state,
    "mode_energies": [v0, v1],
    "mode_couplings": [g0, g1],
    "up_energy": 0.5 * w0,
    "down_energy": -0.5 * w0,
    "energy_cutoff": 10,
    "dt": dt,
}
N = 10000  # number of timesteps
p = pyket.Params()
p.load_dict(d)
h = pyket.H(p)

up_list = np.zeros(N)
down_list = np.zeros(N)
num_states = np.zeros(N)
t = np.zeros(N)
print(h.get_state())
for i in range(N):
    h.grow()
    h.run_step(complex(0, -dt))
    up, down = h.get_spin_pop()
    num_states[i] = h.get_size()
    t[i] = dt * i
    up_list[i], down_list[i] = up, down
fig, ax = plt.subplots(
    2, 1, figsize=(9, 4)
)  # Create a figure containing a single axes.
ax[0].plot(t, up_list, label="up")
ax[0].set_title("up and down vs time")
ax[0].set_ylabel("Probability")
ax[0].plot(t, down_list, label="down")
ax[0].legend()
ax[1].plot(t, num_states, label="number of states")
ax[1].set_ylabel("Number of States")
ax[0].set_xlabel("Time (arbitrary)")
fig.suptitle(f"w0={w0},v0={v0},v1={v1},g0={g0},g1={g1}")
plt.savefig("plot_2_modes." + f"w0={w0}.v0={v0}.v1={v1}.g0={g0}.g1={g1}.png")
