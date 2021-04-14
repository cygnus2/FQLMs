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

def find_level_crossing(s1, s2):
    """ Returns the value of lambda where the level crossing occurs,
        uses interpolation.

        Attention: Assumes a single crossing (which should be reasonable?).
    """
    from scipy.interpolate import interp1d
    from scipy.optimize import minimize_scalar
    f1 = interp1d(s1[0], s1[1])
    f2 = interp1d(s2[0], s2[1])

    eps = 0.0000000
    x = np.array([np.max([s1[0][0], s2[0][0]]), np.min([s1[0][-1], s2[0][-1]])])
    x[0] += eps
    x[1] -= eps
    print(x)

    if np.sign(f1(x[-1])-f2(x[-1])) ==  np.sign(f1(x[0])-f2(x[0])):
        return None, (f1, f2)

    mfun = lambda x: np.abs(f1(x) - f2(x))
    res = minimize_scalar(mfun, args=(), method='bounded', tol=None, bounds=x)
    return res.x, (f1, f2)



datasets = {
    (2,2,2) : {
        '$W = [000]$' : [0,0,0],
        '$W = [001]$' : [0,0,1],
    },
    (2,2,4) : {
        '$W = [000]$' : [0,0,0],
        '$W = [001]$' : [0,0,1]
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
        'gauge_particles' : 'fermions'
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
                lam, gs,
                ls=style[latt]['ls'], marker='.', color='gray', alpha=0.2,
                label=''
            )

        lam_c, (f1, f2) = find_level_crossing((lams[0], energies[0]), (lams[1], energies[1]))
        xf = np.linspace(-0.5, 0.99, 100)
        ax.plot(xf, f1(xf), ls='-')
        ax.plot(xf, f2(xf), ls='-')
        if lam_c is not None:
            ax.axvline(lam_c, ls='--', color='black', lw=0.5)

        # if len(energies) > 1:
        #     # m = energies[0] > energies[1]
        #     # ax.plot(lam[0][~m], chis[0][~m], **style, label='$L=[{:d}{:d}{:d}], W = [000]$'.format(*param['L']))
        #     # ax.plot(lams[m], chis[1][m], **style, label='$L=[{:d}{:d}{:d}], W = [001]$'.format(*param['L']))
        #     ax.plot(lams[0], chis[0], **style[latt], label='$L=[{:d}{:d}{:d}], W=[000]$'.format(*param['L']))
        # else:
        #     ax.plot(lams[0], chis[0], **style, label='$L=[{:d}{:d}{:d}], W=[000]$'.format(*param['L']))


    # ax.set_ylim(-0.1, 2)
    ax.set_xlabel("$\\lambda$")
    ax.set_ylabel("$\\chi_{\\rm F}$")
    ax.legend(loc='upper left')

    fig.tight_layout()
    fig.savefig('plots/mult_fidelity.png')
