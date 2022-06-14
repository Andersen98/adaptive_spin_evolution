# Outline of code structure

Below we outline the `src/hamiltonian/configuration.hpp` and the `src/main.cpp`. 

# Why are we doing this?

Our goal is to create an interface so that we can focus on running simulations and analysis. Having documentaion on just the interface with a few examples will allow us to generally know what we have to do to interface with the code.


## configuration class

We first borrow heavily from the `main.cpp` file, which gives an interface to the configuration class. Here are a few of the class methods.

StateKet (This is a building block class. It represents 1 mode of the field)
- StateKet(::Integer) --> ::StateKet
- get_mode(::Integer) --> 
- set_mode(::Integer)
- get_spin() --> ::Spin (which is just a Bool in our implimentaion)
- set_spin(::Spin) --> None
- set_amp ??




## Hamiltonian Class

First we look at the interface make in `src/hamiltonian/hamiltonian.cpp`. The hamiltonian class is doing a lot and should probably be split up. Below I try to summarize it:

> [!summary]
>
> We can break things down by collecting them: 
> 	- *collections*
> 		- actions: we collect ideas of *processes* that are called in the simulation
> 		- static objects: Things like the geometry and basis sets
> 		- mutable object: Things that *change* as we apply an action
> 		- observables/results: Things like the spectrum or spin density matrix evoluion as a function of time

Now we can just add interfaces/objects/functions to each of the collections.

## mutable objects

What makes this code interesting is that the **active space evolves in time**. This hits on a key motivation for this code.

- Master equations tell us key properties of a system, but only if the full system's states can be split into a reservoior, base system (like a molecule), and an interaction.
- Sometimes it is really hard to know how to split up the interaction from the reservour
- *The resulting active space from our code is desireable because:
  - it tells us how to split our space up for master equations
  - Ex: now one could construct a master equation that well describes the dynamics of a spin system interacting in a non trivial way with the field. 




## Actions


The core of this code is to grow out the active space so that it describes the dynamics of observables like spin density matrix over time. Thus this code's *'purpose*' is to do an action.

**Abstract Actions:*
- time propagate
- grow active space

**Implimentation Actions**

We want to implement the **grow active space** action. We choose to implement it as:
- grow active space = add states (to active space) that have a transition strength above some threshold (where transition strength is proportional to the propagator).
  - partition *the set of states connected by the transition matrix*: partition by if their transition strength is above some threshold
	- make a **mutable** copy of our current (sorted state)
	- propagate some state in our active space: similar to rotating a 2D arrow in a 3D space. The resulting rotation may have taken us out of plane.
		- we **fix** our initial state and only add to the copy
		- to speed up computation, we step through the sorted list from the previous run, as well as a sorted list of couplings. This lets us only propagate above  threshold.	
	- sort the **mutated state** 
		- resulting state into components who's transitio is above some threshold:
	- adding state from that pass the threshold from the sorted result state
		- we start at the strongest connected state (top of list) and step down the list until we find a state below threshold
		- there are two copies of the state: 
			- one that we sort on: 
			- one that we add the 'new' states to : which are just a projection of the propagated state onto the states 'space to our old active space'

### Illustrating with a toy problem


Going back to the 2D arrow example, this is like applying a propagation and (suppose it is in 6D space). The propagation rotates the 2D arrow into the 1st, 2nd, 3rd, 4th, and 5th dimension.  We look at the propagation strength of the 4th dimension is the only one that passed threshold. We add the 4th dim propagation to our initial 2d arrow so the resulting state is in 3D (1st, 2nd, and 4th dimensional space). We then repeat the process until we find a suitable subspace that d
escribes the motion of the arrow. 

## Sifting out the interface

For now I am choosing to interface with the code implimentation like the rotating 2D arrow in 6D space. 

### Group 1

Given [an initial state that spans the active space, an ordering of components by their magnitude, a symmetric transition matrix] I
- I return the part of the propagated vector that passed the threshold (i.e rotated outside the current active space 'in a significant' way)

Some actual methods that help us do this are:
- `bool grow_configuration_space(int idx`
- `void add_connection(state_ket &k, state_vector &psi_mixed)`
  - `ket_pair get_connected_states(const state_ket &k, const int mode)`
  - abstract function `T(s in {set of dims})-> transition strength in reals`
  - a way to make `{ (orderd set of s in dims): where T is monotonically increasing}`




