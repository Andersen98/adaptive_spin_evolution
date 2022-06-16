
# Designing A “Simulation Algebra” for Optical Cavity Simulation Parameters

We have already established that a cavity code usually separate into #concept\s:

- #actions
- #static-objects
- #mutable-objects
- #result-objects

If we come up with some well defined rules on how the #concepts mix and add, we could come up with a sort-of algebra that lets us build complicated simulations from simple building blocks.



## First Attempt

One could imagine that an #action can act as an operator on #objects.

#action * #object = #object 

However the above statement fails if we consider:
> #action:: #transform-state * #result-objects:\: #Field-Component

Some cases where the above _does_ work would be:
> #action:: #transform-state * #mutable-object :\: #State-Vector


## Focus on defining our #action

In the #adaptive_spin_configuration project, the main #action is #get-primary-propagation-components. The action just gets the components of a state that have a transition rate +amplitude above a certain magnitude. 

To use our #action on other libraries, we would just need the transition matrix to be well defined. 




 