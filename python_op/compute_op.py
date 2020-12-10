import numpy as np
import h5py as hdf
import matplotlib.pyplot as plt

import sys
sys.path.append('../python_gauss_lattice/')
from gauss_lattice import HamiltonianBuilder

# Load the stuff for the single states.
ops = np.load('../python_notebooks/counted_plaquettes.npy')
labels = [
    '[xy] $\\langle UU^\\dagger \\rangle$',
    '[xy] $\\langle U^\\dagger U \\rangle$',
    '[yz] $\\langle UU^\\dagger \\rangle$',
    '[yz] $\\langle U^\\dagger U\\rangle$',
    '[xz] $\\langle UU^\\dagger \\rangle$',
    '[xz] $\\langle U^\\dagger U\\rangle$',
]

# Read the GS.
datadir = '../python_data/states_histograms/'
L = [2,2,4]
datafile = datadir+'results_bosons_wx_4-wy_4-wz_2_2x2x4.hdf5'

# Count number of
builder = HamiltonianBuilder({'L':L}, [], silent=True)
plaquettes = builder.get_plaquette_list(separate_lists=True)
n_plaq = [len(p) for p in plaquettes]

# Count number of lambda values.
lambdas = []
with hdf.File(datafile, 'r') as f:
    for g in f:
        if g.startswith("eigenstates"):
            lambdas.append(float(g.split('_')[-1]))
s = np.argsort(lambdas)
lambdas = np.array(lambdas)[s]

# Compute the order-parameter.
n_states = 15
order_parameters = np.zeros(shape=(n_states,len(lambdas),len(labels)))
with hdf.File(datafile, 'r') as f:
    for j, lam in enumerate(lambdas):
        ds = 'eigenstates_lam_{:.6f}'.format(lam)
        for k in range(n_states):
            gs = f[ds][...,k]
            order_parameters[k,j,:] = [np.dot(np.conj(gs)*gs / n_plaq[int(i/2)], ops[:,i]) for i in range(len(labels))]

with plt.style.context('seaborn'):
    nr, nc = 5, 3
    fig, ax = plt.subplots(nr,nc)
    ax = ax.flatten()
    fig.set_size_inches(6*nc, 4*nr)

    for k in range(n_states):
        for j in range(len(labels)):
            shift = 0.05
            ax[k].plot(lambdas+j*shift, order_parameters[k,:,j], marker='o' if j%2 else 's', ls='--', label=labels[j])

    # ax.plot(lambdas, np.sum(order_parameters, axis=1)/2, marker='s', ls='', label='\\langle O_{\\rm f}\\rangle')

        ax[k].set_ylim(0,0.5)
        ax[k].set_xlabel('$\\lambda$')
        ax[k].set_ylabel('$\\langle O \\rangle / N_p$')


    ax[0].legend(loc='upper right', fontsize=6, ncol=2)

fig.tight_layout()
fig.savefig('plots/op_shifted.png')
