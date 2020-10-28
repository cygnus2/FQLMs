""" ----------------------------------------------------------------------------

    gl_aux.py - LR, March 2020

    Holds some auxiliary stuff.

---------------------------------------------------------------------------- """
import time, yaml, sys
import numpy as np
import h5py as hdf
import numpy as np
from itertools import product
import datetime as dt

def timestamp():
    return dt.datetime.now().strftime("[%H:%M:%S]")

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


def bin_str(state, L):
     return ("{:0"+str(L)+"b}").format(state)

def print_2D_state(state, L, spins=False):
    """ Dumps a 2D lattice.
    """
    sx, sy = 2, 1
    blx = 2*(sx+1)

    Lx, Ly = L

    bstr = bin_str(state, Lx*Ly*2)[::-1]

    if spins:
        left, right = "<", ">"
        up, down = "v", "^"
    else:
        left, right = '●', '○'
        up, down = '○', '●'

    bstr = bstr.replace('0', right)
    bstr = bstr.replace('1', left)

    # Loop through rows.
    lines = []
    for r in range(Ly):
        lx = ''
        for c in range(Lx):
            i = 2*r*Lx+c*2
            lx += "-"+"-"*sx + bstr[i] + "-"*sx
        lines.append(lx)

        # ---

        for k in range(sy):
            lines.append(("|"+(" "*blx)[:-1])*Lx)

        ly = ''
        for c in range(Lx):
            j = 2*r*Lx+c*2 + 1
            ly += bstr[j]+(" "*blx)[:-1]
        lines.append(ly.replace(right, up).replace(left, down))

        for k in range(sy):
            lines.append(("|"+(" "*blx)[:-1])*Lx)

    for line in lines[::-1]:
        print(line)


def winding_tag(ws, labels=['x', 'y', 'z'], shift=None):
    """ Returns the naming convention of the winding datasets.

        The shift is a lattice configuration L = [Lx,Ly,Lz], such
        that the HDF5 datasets may be resolved.
    """
    if shift is not None:
        ws = _winding_shift(shift, ws)
    wtag = ''
    for k in range(len(ws)):
        wtag += 'w{:s}_{:d}-'.format(labels[k], ws[k])
    return wtag[:-1]

def _winding_shift(L, ws):
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
    return ws + shift

def read_winding_sector(L, ws, debug=True, basedir='./output/', filename=None):
    """ Takes in a parameter dictionary and reads in the appropriate states
        for the specified winding sector.
    """
    # The winding sectors are labelled differently in the HDF5 file (for
    # convenience reasons) so we have to relabel them here.
    ws_shifted = _winding_shift(L, ws)
    if debug:
        print(f'Supplied winding numbers {ws} are mapped to {tuple(ws_shifted)}')

    # Read and return the appropriate list.
    if not filename:
        filename=basedir+file_tag(L, filetype='hdf5')
    with hdf.File(filename, 'r') as f:
        winding_states = f[winding_tag(ws_shifted)][...]
    return winding_states, winding_tag(ws_shifted)



def read_all_states(L, merged=True, basedir='./output/', filename=None):
    """ Takes in a parameter dictionary and reads in the appropriate states
        for the specified winding sector.
    """
    # Read and return the appropriate list.
    if not filename:
        filename =  basedir+file_tag(L, filetype='hdf5')
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
