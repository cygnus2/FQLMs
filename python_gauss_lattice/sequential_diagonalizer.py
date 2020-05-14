""" ----------------------------------------------------------------------------

    sequential_diagonalizer.py - LR, May 2020

    Run script for diagonalizing the Hamiltonian for the gauss lattice - goes
    through all winding sectors one-by-one.

---------------------------------------------------------------------------- """
from gauss_lattice import HamiltonianBuilder
from gauss_lattice.aux import size_tag, timeit, read_all_states, read_sequential_spectrum, write_simple_spectrum
import numpy as np
import h5py as hdf


@timeit
def hamiltonian_diagonalization(ham, full_diag=False, **kwargs):
    """ This method just exists to be able to time the routine efficiently.
    """
    if full_diag:
        return ham.full_diagonalization()
    return ham.compute_lower_spectrum(**kwargs)


# This sets the parameters fof the calculation (everything else is fixed).
param = {
    'L' : [2,2,4],
    'gauge_particles' : 'bosons'
}
filename = 'output/SEQUENTIAL_spectrum_'+param['gauge_particles']+'_'+size_tag(param['L']) + '.hdf5'
with hdf.File(filename, 'w') as f:
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
    n_eigenvalues = min(100, builder.n_fock//2)
    spectrum = hamiltonian_diagonalization(ham, full_diag=False, n_eigenvalues=n_eigenvalues, which='BE')

    # Save into a dataset in the HDF5 file.
    with hdf.File(filename, 'a') as f:
        f.create_dataset(ws, data=np.array(spectrum))


# Clean up (keep the HDF5 file, but also produce an easier to read spectrum file).
spectrum = read_sequential_spectrum(filename)
write_simple_spectrum(spectrum, filename)
