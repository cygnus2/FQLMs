""" ----------------------------------------------------------------------------

    sequential_diagonalizer.py - LR, May 2020

    Run script for diagonalizing the Hamiltonian for the gauss lattice - goes
    through all winding sectors one-by-one.

---------------------------------------------------------------------------- """
from gauss_lattice import HamiltonianBuilder
from gauss_lattice.aux import param_tag, file_tag, size_tag, timeit, read_all_states, write_simple_spectrum
import numpy as np
import h5py as hdf



def convert_sequential_spectrum(filename, which='BE'):
    """ Takes in a HDF5 file which contains the spectra of several winding
        sectors and combines them in a single list.

        Returns the lowest part of the entire spectrum, up to the value where
        it is sure that it is the lowest value. This cutoff is given by the lowest
        maximum of the sub spectra.

        Disclaimer: only tested for which='BE', but this should work either way.
    """
    spectrum = []
    with hdf.File(filename, 'r') as f:
        min_max = np.inf
        for ds in f:
            sub_spectrum = np.array(sorted(f[ds][...]))
            if which == 'BE':
                cutoff = len(sub_spectrum) // 2
            else:
                cutoff = len(sub_spectrum)
            if cutoff:
                m = np.max(sub_spectrum[:cutoff])
                if m < min_max:
                    min_max = m
                spectrum += sub_spectrum.tolist()

    spectrum = np.array(sorted(spectrum))
    return spectrum[spectrum<min_max]



@timeit
def hamiltonian_diagonalization(ham, **kwargs):
    """ This method just exists to be able to time the routine efficiently.
    """
    return ham.diagonalize(**kwargs)


# This sets the parameters fof the calculation (everything else is fixed).
param = {
    'L' : [2,2,2],
    'J' : -1.0,
    'lambda' : 0,
    'gauge_particles' : 'fermions',

    'ev_type' : 'BE',
    'n_eigenvalues' : 100,

    'store_hamiltonian' : True,
}


spectrum_file = 'output/SEQUENTIAL_spectrum_' + param_tag(param) + '.hdf5'
with hdf.File(spectrum_file, 'w') as f:
    pass

if param['store_hamiltonian']:
    hamiltonian_file = 'output/SEQUENTIAL_hamiltonian_' + size_tag(param['L']) + '.hdf5'
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


    if param['store_hamiltonian']:
        with hdf.File(hamiltonian_file, 'a') as f:
            f.create_dataset(ws, data= np.array([ham.col, ham.row, ham.data]))


    # Diagonalization.
    spectrum = hamiltonian_diagonalization(ham,
        full_diag=False,
        J = param['J'],
        lam = param['lambda'],
        gauge_particles = param['gauge_particles'],
        n_eigenvalues = max(1, min(param['n_eigenvalues'], builder.n_fock//2)),
        which = param['ev_type']
    )

    # Save into a dataset in the HDF5 file.
    with hdf.File(spectrum_file, 'a') as f:
        f.create_dataset(ws, data=np.array(spectrum))

# Clean up (keep the HDF5 file, but also produce an easier to read spectrum file).
spectrum = convert_sequential_spectrum(spectrum_file)
write_simple_spectrum(spectrum, spectrum_file.replace('.hdf5', '.dat'))
