""" ----------------------------------------------------------------------------

    gl_aux.py - LR, March 2020

    Holds some auxiliary stuff.

---------------------------------------------------------------------------- """
import time, yaml, sys
import numpy as np
import h5py as hdf
import numpy as np
from itertools import product


def size_tag(L):
    stag = ''
    for k in range(len(L)):
        stag += '{:d}x'.format(L[k])
    return stag[:-1]

def file_tag(L, filetype='hdf5'):
    """ Returns the naming convention of the winding datasets.
    """
    return 'winding_states_' + size_tag(L) + '.' + filetype

def param_tag(param, precision=2):
    """ Fixed naming convention.
    """
    tag = param['gauge_particles'] + "_" + size_tag(param['L'])

    J = param.get('J')
    if J:
        tag += ('_J{:.'+str(precision)+'f}').format(J)

    lam = param.get('lam')
    if lam:
        tag += ('_lam{:.'+str(precision)+'f}').format(lam)

    return tag


def timeit(logger=None):
    def wrap(method):
        def timed(*args, **kw):
            ts = time.time()
            result = method(*args, **kw)
            te = time.time()
            text = '%r  %2.2f ms' % (method.__name__, (te - ts) * 1000)
            if logger:
                logger.info(text)
            else:
                print(text)
            return result
        return timed
    return wrap


def print_2D_state(state, L):
    """ Dumps a 2D lattice.
    """
    sx, sy = 2, 1
    blx = 2*(sx+1)

    bstr = bin_str(state, L=L**2*2)[::-1]
    bstr = bstr.replace('0', '○')
    bstr = bstr.replace('1', '●')

    # Loop through rows.
    lines = []
    for r in range(L):
        lx = ''
        for c in range(L):
            i = 2*r*L+c*2
            lx += "-"+"-"*sx + bstr[i] + "-"*sx
        lines.append(lx)

        # ---

        for k in range(sy):
            lines.append(("|"+(" "*blx)[:-1])*L)

        ly = ''
        for c in range(L):
            j = 2*r*L+c*2 + 1
            ly += bstr[j]+(" "*blx)[:-1]
        lines.append(ly)

        for k in range(sy):
            lines.append(("|"+(" "*blx)[:-1])*L)

    for line in lines[::-1]:
        print(line)


def winding_tag(ws, labels=['x', 'y', 'z']):
    """ Returns the naming convention of the winding datasets.
    """
    wtag = ''
    for k in range(len(ws)):
        wtag += 'w{:s}_{:d}-'.format(labels[k], ws[k])
    return wtag[:-1]


def read_winding_sector(L, ws, debug=True):
    """ Takes in a parameter dictionary and reads in the appropriate states
        for the specified winding sector.
    """
    # The winding sectors are labelled differently in the HDF5 file (for
    # convenience reasons) so we have to relabel them here.
    if len(L) == 2:
        shift = np.array(L[::-1]) // 2
    elif len(L) == 3:
        shift = np.array([
            L[1]*L[2] // 2,
            L[0]*L[2] // 2,
            L[0]*L[1] // 2
        ])
    else:
        raise NotImplementedError('Dimension not implemented!')
    ws_shifted = ws + shift

    if debug:
        print(f'Supplied winding numbers {ws} are mapped to {tuple(ws_shifted)}')

    # Read and return the appropriate list.
    filename='output/'+file_tag(L, filetype='hdf5')
    with hdf.File(filename, 'r') as f:
        winding_states = f[winding_tag(ws_shifted)][...]
    return winding_states



def read_all_states(L, merged=True):
    """ Takes in a parameter dictionary and reads in the appropriate states
        for the specified winding sector.
    """
    # Read and return the appropriate list.
    filename='output/'+file_tag(L, filetype='hdf5')
    states = []
    with hdf.File(filename, 'r') as f:
        for ws in f:
            if merged:
                states += list(f[ws][...])
            else:
                states.append([ws, list(f[ws][...])])
    return states


def read_sequential_spectrum(filename):
    spectrum = []
    with hdf.File(filename, 'r') as f:
        for ds in f:
            spectrum += f[ds][...].tolist()
    return sorted(spectrum)


def write_simple_spectrum(spectrum, filename):
    with open(filename, 'w') as f:
        for e in spectrum:
            f.write('{:.8f},{:.8f}\n'.format(e.real, e.imag))




def winding_sectors(L, tag=False):
    """ A generator for all winding number sectors.
    """
    if len(L) == 2:
        winding_numbers = tuple(np.array(L[::-1])+1)
    elif len(L) == 3:
        winding_numbers = (
            L[1]*L[2] + 1,
            L[0]*L[2] + 1,
            L[0]*L[1] + 1
        )
    else:
        raise NotimplementedError('Only 2D and 3D lattices are allowed.')

    sectors = product(*map(lambda n: range(n), winding_numbers))
    for sector in sectors:
        if tag:
            yield winding_tag(sector)
        else:
            yield sector


def load_config(filename):
    """ Reads a config from a yaml file, with appropriate chekcs.
    """
    if filename is not None:
        # Loads parameters from file.
        with open(filename, 'r') as f:
            try:
                return yaml.safe_load(f)
            except yaml.YAMLError as exc:
                print(exc)
                raise yaml.YAMLError()
    else:
        sys.exit('fatal: no input file specified')
