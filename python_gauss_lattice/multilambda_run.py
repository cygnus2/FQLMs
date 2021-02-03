""" ----------------------------------------------------------------------------

    multilambda_runs.py - LR, October 2020



---------------------------------------------------------------------------- """
from gl_simulation import GLSimulation
from gauss_lattice import ParallelHamiltonianBuilder, HamiltonianBuilder
import numpy as np


sim = GLSimulation()

ham = sim.read_hamiltonian(ham_name=sim.ws)
if not ham:
    # Produce states.
    states = sim.read_states()
    if not len(states):
        states = sim.find_states()

    ham = sim.construct_hamiltonian(states, builder_type=ParallelHamiltonianBuilder)
    for k in range(len(ham.data)):
        print("({:d},{:d}) // {:.2f}".format(ham.row[k]+1, ham.col[k]+1, ham.data[k]))

    if sim.store_ham:
        sim.store_hamiltonian(ham, label=sim.ws, attrs={'n_fock':len(states)})


lambdas = np.linspace(*sim.param['lambdas']) if 'lambdas' in sim.param else [sim.param['lambda']]
sim.run_lambda_loop(lambdas, ham)
sim.push_notification('Done with multiparameter run.')
