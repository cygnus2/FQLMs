""" ----------------------------------------------------------------------------

    gs_fidelity.py - LR, March 2021

    Here we want to explore the GS fidelity for 3D QLM systems.

---------------------------------------------------------------------------- """
import matplotlib.pyplot as plt
import numpy as np
from copy import copy


import sys
sys.path.append('../python_data/')
from ed_result_wrapper import EDResult

# Read data.
param = {
    'L' : [2,2,2],
    'gauge_particles' : 'bosons',
}

datafiles = [None, 'results_{:s}_all-ws_2x2x2_fine.hdf5'.format(param['gauge_particles'])]
results = [EDResult(param, datadir="data/", datafile=df) for df in datafiles]

# Compute gs fidelity as function of lambda.
fidelities = []
chis = []
entropies, iprs = [], []
for result in results:
    lambdas = np.sort(np.array(list(set(result.ev['lambda']))))
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

    fidelities.append([lambdas[1:], fidelity])
    chis.append([lambdas[1:], chi])


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
    fig, ax = plt.subplots(4,1)
    fig.set_size_inches(5,12)

    bx = ax[2].twinx()

    # Skip RK point.
    c = -1#-35

    for k, (x, y) in enumerate(fidelities):
        ax[0].plot(x[:c], y[:c], label='$d\\lambda = {:.5f}$'.format(x[1]-x[0]))

    for k, (x, y) in enumerate(chis):
        ax[1].plot(x[:c], y[:c], label='$d\\lambda = {:.5f}$'.format(x[1]-x[0]))

    for k, (x, y) in enumerate(entropies):
        ax[2].plot(x[:c], y[:c], label='$d\\lambda = {:.5f}$'.format(x[1]-x[0]))
        bx.plot(x[1:c], np.diff(y[:c])/np.diff(x[:c]), ls='--', label='dS/dlambda')

    for k, (x, y) in enumerate(iprs):
        ax[3].plot(x[:c], y[:c], label='$d\\lambda = {:.5f}$'.format(x[1]-x[0]))


    bx.axis('off')

    ax[0].set_xlabel('$\\lambda$')
    ax[0].set_ylabel('$|\\langle\\psi_0(\lambda)|\\psi_0(\lambda+d\\lambda)\\rangle|$')
    ax[0].legend(loc='lower left')
    ax[0].set_title('GS fidelity [{:s}, {:d}x{:d}x{:d}]'.format(param['gauge_particles'], *param['L']))

    ax[1].set_xlabel('$\\lambda$')
    ax[1].set_ylabel('$\\chi_F$')
    ax[1].legend(loc='upper left')
    ax[1].set_title('fidelity susceptibility')

    ax[2].set_xlabel('$\\lambda$')
    ax[2].set_ylabel('$S_1$')
    ax[2].legend(loc='upper left')
    ax[2].set_title('Shannon entropy (+ numerical derivative)')

    ax[3].set_xlabel('$\\lambda$')
    ax[3].set_ylabel('IPR $ = -\\ln(S_2)$')
    ax[3].legend(loc='upper right')
    ax[3].set_title('inverse participation ration (IPR)')


    for a in ax:
        a.set_xlim(-3,1)

    fig.tight_layout()
    fig.savefig('gs_fidelity_{:s}_{:d}x{:d}x{:d}.pdf'.format(param['gauge_particles'], *param['L']))
