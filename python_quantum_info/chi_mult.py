""" ----------------------------------------------------------------------------

    gs_fidelity.py - LR, March 2021

    Here we want to explore the GS fidelity for 3D QLM systems.

---------------------------------------------------------------------------- """
import matplotlib.pyplot as plt
import numpy as np
from copy import copy
from aux import get_result


def get_overlap_chi_from_result(result):
    """ Computes the fidelity susceptibilty via overlap.
    """
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

    gs = result.ev[result.ev['N']==0].sort_values(by='lambda')['spectrum'].values
    return np.array(lambdas), chi, np.array(gs)


datasets = {
    (2,2,2) : {
        '$W = [000]$' : [0,0,0],
        # '$W = [001]$' : [0,0,1],
    },
    (2,2,4) : {
        '$W = [000]$' : [0,0,0],
        # '$W = [001]$' : [0,0,1]
    }
}
style = {
    (2,2,2): {'marker':'o', 'ls':'-'},
    (2,2,4): {'marker':'s', 'ls':'--'},
}

with plt.style.context('seaborn'):
    fig, ax = plt.subplots()

    # Read data.
    param = {
        'gauge_particles' : 'bosons'
    }

    c = -1
    for latt, winding_sectors in datasets.items():
        param['L'] = latt
 
        lams, chis, energies = [], [], []
        for label, ws in winding_sectors.items():
            param['winding_sector'] = ws
            result, _ = get_result(param, datadir='data/', use_pickle=(True,True))

            lam, chi, gs = get_overlap_chi_from_result(result)
            print(len(lam), " // ", len(chi), " // ", len(gs))
            lams.append(lam)
            energies.append(gs)
            chis.append(chi)
            ax.plot(
                lam[:-1], chi,
                ls=style[latt]['ls'], marker='o',
                label='L [{:d}{:d}{:d}]'.format(*latt)
            )

        # ax.plot(lams[0][:-1], chis[0], **style, label='$L=[{:d}{:d}{:d}], W=[000]$'.format(*param['L']))


    ax.set_title(param['gauge_particles'])
    ax.set_ylim(-0.1, 2)
    ax.set_xlabel("$\\lambda$")
    ax.set_ylabel("$\\chi_{\\rm F}$")
    ax.legend(loc='upper left')

    fig.tight_layout()
    fig.savefig('plots/mult_fidelity.png')
