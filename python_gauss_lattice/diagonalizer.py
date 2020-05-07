""" ----------------------------------------------------------------------------

    diagonalizer.py - LR, May 2020

    Run script for diagonalizing the Hamiltonian for the gauss lattice.

---------------------------------------------------------------------------- """
from gauss_lattice import HamiltonianBuilder
from gauss_lattice.aux import file_tag, winding_tag, timeit
import numpy as np
import h5py as hdf


def read_winding_sector(param, filename):
    """ Takes in a parameter dictionary and reads in the appropriate states
        for the specified winding sector.
    """
    ws = param.get('winding_sector')
    if not ws:
        raise KeyError("No winding sector specified!")

    # The winding sectors are labelled differently in the HDF5 file (for
    # convenience reasons) so we have to relabel them here.
    L = param['L']
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
    print(ws)
    print(ws_shifted)
    with hdf.File(filename, 'r') as f:
        winding_states = f[winding_tag(ws_shifted)][...]
    return winding_states


@timeit
def hamiltonian_construction(builder):
    """ This method just exists to be able to time the routine efficiently.
    """
    return builder.construct()

@timeit
def hamiltonian_diagonalization(ham, *args, **kwargs):
    """ This method just exists to be able to time the routine efficiently.
    """
    return ham.compute_lower_spectrum(*args, **kwargs)


# This sets the parameters fof the calculation (everything else is fixed).
param = {
    'L' : [2,2,2],
    'winding_sector' : (0,0,0),
    'J' : 1.0
}

winding_states = read_winding_sector(param, filename='output/'+file_tag(param['L'], filetype='hdf5'))

# Set up the builder object.
builder = HamiltonianBuilder(
    param,
    states=winding_states
)

# Get the hamiltonian.
ham = hamiltonian_construction(builder)

# Diagonalization.
eigenvalues = hamiltonian_diagonalization(ham, n_eigenvalues=50, which='BE', dense=False)

# Some I/O.
print(eigenvalues)
ham.store_results(filename='output/spectrum_'+file_tag(param['L'], filetype='dat'))
