""" ----------------------------------------------------------------------------

    h_constructor.py - LR, October 2020

    Constructs and stores a Hamiltonian.

---------------------------------------------------------------------------- """
from gauss_lattice import LowEnergyHamiltonianBuilder, ParallelHamiltonianBuilder, HamiltonianBuilder
from gl_simulation import GLSimulation


sim = GLSimulation()

states = sim.read_states()
if not states:
    raise ValueError('Fatal: no states found.')

# Here we decide which object we used to construct the Hamiltonian!
ham = sim.construct_hamiltonian(states, builder_type=ParallelHamiltonianBuilder)
if sim.store_ham:
    sim.store_hamiltonian(ham, label=sim.ws, attrs={'n_fock':len(states)})

sim.push_notification('Done with Hamiltonian construction.')
