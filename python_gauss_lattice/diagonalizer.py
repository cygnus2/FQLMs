""" ----------------------------------------------------------------------------

    diagonalizer.py - LR, May 2020

    Run script for diagonalizing the Hamiltonian for the gauss lattice.

---------------------------------------------------------------------------- """
from hamiltonian_builder import HamiltonianBuilder
from gl_aux import file_tag, winding_tag, timeit
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
    ws_shifted = np.array(ws) + np.array(param['L'])
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
    'L' : [2, 2, 2],
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
eigenvalues = hamiltonian_diagonalization(ham)

# Some I/O.
print(eigenvalues)
ham.store_results(filename='output/lower_spectrum_2x2x2.dat')
