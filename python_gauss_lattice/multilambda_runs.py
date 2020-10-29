""" ----------------------------------------------------------------------------

    multilambda_runs.py - LR, October 2020



---------------------------------------------------------------------------- """
from simulation import GLSimulation
import numpy as np


sim = GLSimulation()

# Produce states.
states = sim.read_states()
if not states:
    raise ValueError('Fatal: no states found.')

ham = sim.read_hamiltonian(ham_name=sim.ws)
if not ham:
    ham = sim.construct_hamiltonian(states)
    if sim.store_ham:
        sim.store_hamiltonian(ham, label=sim.ws, attrs={'n_fock':len(states)})


lambdas = np.linspace(*sim.param['lambdas'])
for i, l in enumerate(lambdas):
    sim.log('[{:d} / {:d}] diagonalizing Hamiltonian for lambda={:.4f}'.format(i+1, len(lambdas), l))

    # Diagonalization.
    results = sim.diagonalize_hamiltonian(ham, lam=l)
    attrs = {
        'lambda' : l
    }
    if not sim.compute_eigenstates:
        sim.store_data(results, ds_name='spectrum_lam_{:6f}'.format(l), attrs=attrs)
    else:
        sim.store_data(results[0], ds_name='spectrum_lam_{:6f}'.format(l), attrs=attrs)
        sim.store_data(results[1], ds_name='eigenstates_lam_{:6f}'.format(l), attrs=attrs)


sim.push_notification('Done with multiparameter run.')
