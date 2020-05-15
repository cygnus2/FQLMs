""" ----------------------------------------------------------------------------

    sequential_diagonalizer.py - LR, May 2020

    Run script for diagonalizing the Hamiltonian for the gauss lattice - goes
    through all winding sectors one-by-one.

---------------------------------------------------------------------------- """
from gauss_lattice import HamiltonianBuilder
from gauss_lattice.aux import param_tag, file_tag, timeit, read_all_states, read_sequential_spectrum, write_simple_spectrum
import numpy as np
import h5py as hdf


@timeit
def hamiltonian_diagonalization(ham, **kwargs):
    """ This method just exists to be able to time the routine efficiently.
    """
    return ham.diagonalize(**kwargs)


# This sets the parameters fof the calculation (everything else is fixed).
param = {
    'L' : [2,2,2],
    'J' : -1.0,
    'lambda' : -0,
    'gauge_particles' : 'bosons',

    'ev_type' : 'BE',
    'n_eigenvalues' : 100,

    'store_hamiltonian' : True,
}


spectrum_file = 'output/SEQUENTIAL_spectrum_' + param_tag(param) + '.hdf5'
with hdf.File(spectrum_file, 'w') as f:
    pass

if param['store_hamiltonian']:
    hamiltonian_file = 'output/SEQUENTIAL_hamiltonian_' + param_tag(param) + '.hdf5'
    with hdf.File(hamiltonian_file, 'w') as f:
        pass

# Read all winding sectors from file.
all_winding_sectors = read_all_states(param['L'], merged=False)

total_progress = 0
total_states = sum(list(map(lambda x: len(x[1]), all_winding_sectors)))
for i, winding_sector in enumerate(all_winding_sectors):
    ws, states = winding_sector

    total_progress += len(states)
    print('Diagonalizing sector {:d} of {:d} [{:d} of {:d} states]'.format(i+1, len(all_winding_sectors), total_progress, total_states))

    # Set up the builder object & construct the Hamiltonian.
    builder = HamiltonianBuilder(param, states=states)
    ham = builder.construct()

    # Diagonalization.
    n_eigenvalues = max(1, min(param['n_eigenvalues'], builder.n_fock//2))
    spectrum = hamiltonian_diagonalization(ham, full_diag=False, n_eigenvalues=n_eigenvalues, which=param['ev_type'])

    # Save into a dataset in the HDF5 file.
    with hdf.File(spectrum_file, 'a') as f:
        f.create_dataset(ws, data=np.array(spectrum))

    if param['store_hamiltonian']:
        with hdf.File(hamiltonian_file, 'a') as f:
            f.create_dataset(ws, data=np.array(spectrum))


# Clean up (keep the HDF5 file, but also produce an easier to read spectrum file).
spectrum = read_sequential_spectrum(spectrum_file)
write_simple_spectrum(spectrum, spectrum_file)
