""" ----------------------------------------------------------------------------

    data_io.py - LR, October 2020

    Holds some useful routines to efficiently read data.

---------------------------------------------------------------------------- """
import numpy as np
import pandas as pd
import h5py as hdf


def read_spectrum(datafile):
    return pd.read_csv(datafile, header=None)[0].values


def retrieve_spectrum_pandas(filename):
    """ Returns a pandas dataframe with all speparate spectra.
    """
    df = pd.DataFrame()
    with hdf.File(filename, 'r') as f:
        for ds in f:
            data = {
                'lambda' : float(ds.split('_')[-1]),
                'spectrum' : [f[ds][:]]
            }
            df = df.append(pd.DataFrame(data), ignore_index=True)

    df = df.sort_values(by='lambda')
    return df


def retrieve_spectrum_numpy(filename, cutoff=20, skip=[]):
    """ Returns a 2D numpy array with all the information.
    """
    data = {}
    lambdas = []
    with hdf.File(filename, 'r') as f:
        for ds in f:
            if not any([tag in ds for tag in skip]):
                lam = ds.split('_')[-1]
                data[lam] = f[ds][:]
                lambdas.append(float(lam))
    lambdas = sorted(lambdas)

    spectrum = np.zeros(shape=(len(lambdas),cutoff), dtype=np.complex) * np.nan
    for k, key in enumerate(lambdas):
        d = data['{:.6f}'.format(key)][:cutoff]
        spectrum[k,:len(d)] = d 
    return np.array(lambdas), spectrum


# --------------------------------
# Low energy stuff.

def retrieve_le_spectrum_numpy(filename, grp=None, cutoff=20):
    """ Returns a 2D numpy array with all the information.
    """
    data = {}
    lambdas = []
    with hdf.File(filename, 'r') as f:
        grp = f[grp if grp is not None else '']
        for ds in grp:
            lam = ds.split('_')[-1]
            data[lam] = grp[ds][...]

            lambdas.append(float(lam))
    lambdas = sorted(lambdas)

    spectrum = np.zeros(shape=(len(lambdas),cutoff))
    for k, key in enumerate(lambdas):
        spectrum[k,:] = data['{:.6f}'.format(key)][:cutoff]
    return np.array(lambdas), spectrum


# def retrieve_le_spectrum_numpy(filename, max_level, cutoff=20):
#     """ Returns a 2D numpy array with all the information.
#     """
#     data = {}
#     lambdas = []
#     with hdf.File(filename, 'r') as f:
#         grp = f['ex_{:d}'.format(max_level)]
#         for ds in grp:
#             lam = ds.split('_')[-1]
#             data[lam] = grp[ds][...]

#             lambdas.append(float(lam))
#     lambdas = sorted(lambdas)

#     spectrum = np.zeros(shape=(len(lambdas),cutoff))
#     for k, key in enumerate(lambdas):
#         spectrum[k,:] = data['{:.6f}'.format(key)][:cutoff]
#     return np.array(lambdas), spectrum




def retrieve_le_spectrum_pandas(filename):
    """ Returns a pandas dataframe with all speparate spectra.
    """
    df = pd.DataFrame()
    with hdf.File(filename, 'r') as f:
        for g in f:
            grp = f[g]
            for ds in grp:
                data = {
                    'ex' : int(g.split('_')[-1]),
                    'lambda' : float(ds.split('_')[-1]),
                    'spectrum' : [grp[ds][:]]
                }
                df = df.append(pd.DataFrame(data), ignore_index=True)

    df = df.sort_values(by='lambda')
    return df
