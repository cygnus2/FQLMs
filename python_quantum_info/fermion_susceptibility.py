""" ----------------------------------------------------------------------------

    fermion_susceptibility.py - LR, April 2021

    To check the transition between the winding states.

---------------------------------------------------------------------------- """
import matplotlib.pyplot as plt
import numpy as np
from aux import get_result


def get_chi_from_result(result):
    """ Returns lambdas and chi for a given result.
    """
    ev = result.get_eigenvalues('susceptibility')

    ev = ev[ev['lambda']<1.0] # Sikp RK point.
    lams = np.sort(list(set(ev['lambda'])))

    # Sum over all contributions.
    lcurve = np.zeros_like(lams)
    gs = np.zeros_like(lams)
    for i, lam in enumerate(lams):
        sub = ev[ev['lambda']==lam].sort_values(by='N')
        lcurve[i] = sum(sub['susceptibility'])
        gs[i] = sub[sub['N']==0].iloc[0]['spectrum']
    return lams, lcurve, gs


# Read data.
param = {
    'L' : [2,2,2],
    'gauge_particles' : 'fermions',
}

winding_sectors = {
    '$W = [000]$' : [0,0,0],
    '$W = [001]$' : [0,0,1],
}

with plt.style.context('seaborn-notebook'):
    fig, ax = plt.subplots()

    chis, energies = [], []
    for label, ws in winding_sectors.items():
        param['winding_sector'] = ws
        result, _ = get_result(param, datadir='data/')

        lams, chi, gs = get_chi_from_result(result)
        energies.append(gs)
        chis.append(chi)
        ax.plot(
            lams, chi,
            ls='-', marker='.', color='gray', alpha=0.2,
            label=''
        )

    style = {'ls':'', 'marker':'o'}
    m = energies[0] > energies[1]
    ax.plot(lams[~m], chis[0][~m], **style, label='$W = [000]$')
    ax.plot(lams[m], chis[1][m], **style, label='$W = [001]$')


    ax.set_xlabel("$\\lambda$")
    ax.set_ylabel("$\\chi_{\\rm F}$")
    ax.legend(loc='upper left')

    fig.tight_layout()
    fig.savefig('plots/fermion_fideltiy_222.png')
