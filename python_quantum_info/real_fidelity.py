""" ----------------------------------------------------------------------------

    real_fideliy.py - LR, April 2021

    Here we want to explore the GS fidelity for 3D QLM systems. The "real" refers
    to the computation done with the full spectrum, as opposed to the overlap
    between to neighboring states.

    I checked that those two things do give the same results, however, the former
    let's us investigate also the contributions from different states which could
    give us extra insights?

---------------------------------------------------------------------------- """
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import numpy as np
from copy import copy
import pickle
import os

from aux import get_result

# Read data.
param = {
    'L' : [2,2,2],
    'gauge_particles' : 'fermions',
    'winding_sector' : [1,0,0]
}

# Load result.
result, tag = get_result(param, datadir='data/')

ev = result.get_eigenvalues('susceptibility')
ev = ev[ev['lambda']<1.0] # Sikp RK point.
lams = np.sort(list(set(ev['lambda'])))
lcurve = np.zeros_like(lams)

with plt.style.context('seaborn-notebook'):
    fig, ax = plt.subplots(3,1)
    fig.set_size_inches(6,12)

    lam_field = np.zeros(shape=(len(lams), len(set(ev['N']))))
    for i, lam in enumerate(lams):
        sub = ev[ev['lambda']==lam].sort_values(by='N')
        color = cm.jet(i/len(lams))
        lam_field[i,:] = np.log10(sub['susceptibility'])

        ax[0].plot(
            sub['N'], sub['susceptibility'],
            color=color, ls='', marker='.',
            label='$\\lambda = {:.2f}$'.format(lam)
        )

        s = sum(sub['susceptibility'])
        ax[2].plot(lam, s, marker='o', color=color)
        lcurve[i] = s

    cmap = cm.BuPu
    im = ax[1].imshow(
        lam_field[:,:], origin='lower', aspect='auto', cmap=cmap,
        interpolation='none',
        extent=[0, lam_field.shape[-1], np.min(lams), np.max(lams)]
    )
    ax_divider = make_axes_locatable(ax[1])
    cax = ax_divider.append_axes("right", size="3%", pad="2%")
    fig.colorbar(im, cax=cax, cmap=cmap, orientation='vertical')

    # ax[2].plot(lams, lcurve, color='gray')

    # ---
    ax[0].set_title(tag[1:])

    ax[0].set_yscale('log')
    ax[0].set_xlabel("# of state in basis")
    ax[0].set_ylabel("contribution to fidelity")

    ax[1].set_xlabel("# of state in basis")
    ax[1].set_ylabel("$\\lambda$")

    ax[2].set_xlabel("$\\lambda$")
    ax[2].set_ylabel("$\\chi_{\\rm F}$")

    fig.tight_layout()
    fig.savefig('plots/real_fidelity'+tag+'.png')

exp = {'lambdas':lams, 'chi_f': lcurve}
pickle.dump(exp, open("fidelity"+tag+".pickle", 'wb'))
