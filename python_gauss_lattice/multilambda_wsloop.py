""" ----------------------------------------------------------------------------

    multilambda_wsloop.py - LR, November 2020

    Run script that cycles over all winding-sectors. Makes a bit of a mess
    because there's one file per sector.

    TODO:
     - make it possible to choose winding sectors.

---------------------------------------------------------------------------- """
from gl_simulation import GLSimulation
from gauss_lattice import ParallelHamiltonianBuilder, HamiltonianBuilder
import numpy as np




sim = GLSimulation()
for ws in GLSimulation.winding_sectors(sim.param['L']):
    sim.set_winding_sector(ws, shift=True)

    sim.log('-------------------------------------------')
    sim.log(f'Running ws {ws}')

    ham = sim.read_hamiltonian(ham_name=sim.ws)
    if not ham:

        # Produce states.
        states = sim.read_states()
        if not len(states):
            states = sim.find_states()

        ham = sim.construct_hamiltonian(states, builder_type=ParallelHamiltonianBuilder)
        if sim.store_ham:
            sim.store_hamiltonian(ham, label=sim.ws, attrs={'n_fock':len(states)})

    lambdas = np.linspace(*sim.param['lambdas']) if 'lambdas' in sim.param else [sim.param['lambda']]
    sim.run_lambda_loop(lambdas, ham)
    sim.log('-------------------------------------------')
