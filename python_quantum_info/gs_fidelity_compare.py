""" ----------------------------------------------------------------------------

    gs_fidelity.py - LR, March 2021

    Here we want to explore the GS fidelity for 3D QLM systems.

---------------------------------------------------------------------------- """
import matplotlib.pyplot as plt
import numpy as np
from copy import copy
import pickle


import sys
sys.path.append('../python_data/')
from ed_result_wrapper import EDResult

# Read data.
param = {
    'L' : [2,2,2],
    'gauge_particles' : 'bosons',
}

results = [EDResult(param, datadir="data/", datafile=df) for df in [None]]

# Compute gs fidelity as function of lambda.
fidelities = []
chis = []
entropies, iprs = [], []
for result in results:
    lambdas = np.sort(np.array(list(set(result.ev['lambda']))))
    print(len(lambdas))
    dl = np.diff(lambdas)
    fidelity = np.zeros(shape=(len(lambdas)-1))
    chi = np.zeros_like(fidelity)

    state_left = result.get_eigenstates(lambdas[0])[0,:]
    for k in range(1,len(lambdas)):
        state_right = result.get_eigenstates(lambdas[k])[0,:]
        fidelity[k-1] = np.abs(np.dot(np.conj(state_left), state_right))

        # TODO: Proper normalization required.
        # (some volume factor, probably - this will be interesting for finite size scaling)
        chi[k-1] = (1-fidelity[k-1])*2/dl[k-1]**2
        state_left = copy(state_right)

    fidelities.append([lambdas[:-1], fidelity])
    chis.append([lambdas[:-1]+dl/2, chi])


    # Also compute entanglement.
    ent, ipr = np.zeros(shape=(len(lambdas))), np.zeros(shape=(len(lambdas)))
    for k, lam in enumerate(lambdas):
        gs = result.get_eigenstates(lam)[0,:]
        gs2 = np.conj(gs)*gs
        ent[k] = -np.sum(gs2*np.log(gs2))
        ipr[k] = np.sum(gs2**2)
    entropies.append([lambdas, ent])
    iprs.append([lambdas, ipr])

with plt.style.context('seaborn'):
    fig, ax = plt.subplots()

    # Skip RK point.
    c = -1
    for k, (x, y) in enumerate(chis):
        ax.plot(
            x[:c], y[:c],
            label='approx, $\\epsilon = {:.4f}$'.format(x[1]-x[0])
        )

    # Compare to full calculation.
    pfile = 'fidelity_{:s}_{:d}x{:d}x{:d}.pickle'.format(param['gauge_particles'], *param['L'])
    pdat = pickle.load(open(pfile, 'rb'))
    ax.plot(pdat['lambdas'], pdat['chi_f'], ls='', marker='o', label='full')

    ax.set_xlabel('$\\lambda$')
    ax.set_ylabel('$\\chi_{\\rm F}$')
    ax.legend(loc='upper left')

    fig.tight_layout()
    fig.savefig('gs_fidelity_comparison_{:s}_{:d}x{:d}x{:d}.png'.format(param['gauge_particles'], *param['L']))
