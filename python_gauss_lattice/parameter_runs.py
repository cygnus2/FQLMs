""" ----------------------------------------------------------------------------

    parameter_runs.py - LR, May 2020

    Scans through parameter space starting from a constructed Hamiltonian.

---------------------------------------------------------------------------- """
from gauss_lattice import GaussLatticeHamiltonian, HamiltonianBuilder, GaussLattice
from gauss_lattice.aux import size_tag, timeit, read_all_states, file_tag, load_config
import numpy as np
import h5py as hdf
import argparse, logging

parser = argparse.ArgumentParser(description="Python gauss lattice diagonalizer.")
parser.add_argument('-i', metavar='', type=str, default=None, help='YAML style input file.')
parser.add_argument('-notify', action='store_true', help='Toggles whether a notification should be sent when done. Should only be used by Lukas (sorry).')
args = parser.parse_args()

# This sets the parameters fof the calculation (everything else is fixed).
param = load_config(args.i)

# Set up a logger with a handler for the terminal output.
logger = logging.getLogger('parameter run logger')
logger.addHandler(logging.StreamHandler())
logger.addHandler(logging.FileHandler(param['logfile'], mode='w'))
logger.setLevel(logging.DEBUG)


@timeit(logger=logger)
def hamiltonian_diagonalization(ham, **kwargs):
    """ This method just exists to be able to time the routine efficiently.
    """
    return ham.diagonalize(**kwargs)

# ---
# Retrieves the Hamiltonian and if necessary, constructs it.

input_file = 'output/hamiltonian_' + size_tag(param['L']) + '.npz'
try:
    ham = GaussLatticeHamiltonian.from_scipy_dump(input_file)
    logger.info('Read Hamiltonian from file.')
except FileNotFoundError:
    logger.info('Could not find stored Hamiltonian, attempting to read Fock states.')
    try:
        states = read_all_states(param['L'])
        logger.info('Read Fock states from file, constructing Hamiltonian.')
    except FileNotFoundError:
        logger.info('Could not find stored states, constructing Fock state list from scratch.')
        glatt = GaussLattice(param['L'], state_file=file_tag(L, filetype='hdf5'))
        glatt.find_states()
        states = read_all_states(param['L'])
        logger.info('Constructed states.')
    builder = HamiltonianBuilder(param, states, logger=logger)
    ham = builder.construct()
    ham.store_hamiltonian(input_file)
    logger.info('Created and stored Hamiltonian matrix.')

# ----

spectra = {}
lambdas = np.linspace(*param['lambdas'])
for i, l in enumerate(lambdas):
    logger.info('[{:d} / {:d}] diagonalizing Hamiltonian for lambda={:.4f}'.format(i+1, len(lambdas), l))
    # Diagonalization.
    spectra[l] = hamiltonian_diagonalization(ham,
        full_diag = param.get('full_diag'),
        J = param['J'],
        lam = l,
        gauge_particles = param['gauge_particles'],
        n_eigenvalues = max(1, min(param['n_eigenvalues'], ham.n_fock//2)),
        which = param['ev_type']
    )

    with hdf.File('output/multi_spectrum_'+size_tag(param['L']) + '.hdf5', 'a' if i else 'w') as f:
        f.attrs['lambdas'] = lambdas
        ds = f.create_dataset('lam_{:6f}'.format(l), data=spectra[l])
        ds.attrs['lambda'] = l


if args.notify:
    from gauss_lattice.lr_notify import push_message
    import os, subprocess
    host = subprocess.check_output(['hostname']).strip().decode('UTF-8'),
    push_message(f'Multi-parameter run  @{host[0]} is done!')
