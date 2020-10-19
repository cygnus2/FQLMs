""" ----------------------------------------------------------------------------

    parameter_runs.py - LR, May 2020

    Scans through parameter space starting from a constructed Hamiltonian.

---------------------------------------------------------------------------- """
from gauss_lattice import GaussLatticeHamiltonian, HamiltonianBuilder, GaussLattice
from gauss_lattice.aux import size_tag, timeit, read_all_states, file_tag, load_config, read_winding_sector
import numpy as np
import h5py as hdf
import argparse, logging, os

parser = argparse.ArgumentParser(description="Python gauss lattice diagonalizer.")
parser.add_argument('-i', metavar='', type=str, default=None, help='YAML style input file.')
parser.add_argument('-readonly', action='store_true', help='If this is set, the script only runs if the necessary files are found in the working directory (good for cluster usage if large files would be produced).')
parser.add_argument('-notify', action='store_true', help='Toggles whether a notification should be sent when done. Should only be used by Lukas (sorry).')
args = parser.parse_args()

# This sets the parameters fof the calculation (everything else is fixed).
param = load_config(args.i)
hamiltonian_file = param['working_directory'] + '/SEQUENTIAL_hamiltonian_' + size_tag(param['L']) + '.hdf5'

# Set up a logger with a handler for the terminal output.
logger = logging.getLogger('parameter run logger')
logger.addHandler(logging.StreamHandler())
logger.addHandler(logging.FileHandler(param['working_directory'] + '/' + param['logfile'], mode='w'))
logger.setLevel(logging.DEBUG)


@timeit(logger=logger)
def hamiltonian_construction(builder, *args, **kwargs):
""" This method just exists to be able to time the routine efficiently.
"""
return builder.construct(*args, **kwargs)

@timeit(logger=logger)
def hamiltonian_diagonalization(ham, **kwargs):
    """ This method just exists to be able to time the routine efficiently.
    """
    return ham.diagonalize(**kwargs)


# Check if states exist, if not, construct them.
state_file = param['working_directory']+'/'+file_tag(param['L'], filetype='hdf5')
if not os.path.isfile(state_file):
    logger.info('Could not find stored states, constructing Fock state list from scratch.')
    glatt = GaussLattice(param['L'], state_file=file_tag(param['L'], filetype='hdf5'), basedir=param['working_directory'])
    glatt.find_states()
    logger.info('Constructed states.')
else:
    logger.info('Reading Fock states from file.')

# Read the states.
if param.get('winding_sector'):
    states, ws = read_winding_sector(param['L'], param['winding_sector'], basedir=param['working_directory'])
else:
    states = read_all_states(param['L'], basedir=param['working_directory'])
    ws = None
# ---
# Retrieves the Hamiltonian and if necessary, constructs it.
try:
    with hdf.File(hamiltonian_file, 'r') as f:
        mat = f[ws][...]
    if mat.shape[-1]:
        ham = GaussLatticeHamiltonian(mat[2,:], mat[1,:], mat[0,:], len(states))
    else:
        ham = GaussLatticeHamiltonian([], [], [], 1)
    logger.info('Read Hamiltonian from provided file.')

except (KeyError, OSError):
    # Set up the builder object & construct the Hamiltonian.
    builder = HamiltonianBuilder(param, states=states, logger=logger)
    ham = hamiltonian_construction(builder, param.get('n_threads', 1))

    # If specified, store the Hamiltonian for later use.
    if param['store_hamiltonian']:
        with hdf.File(hamiltonian_file, 'a' if i else 'w') as f:
            ds = f.create_dataset(ws, data= np.array([ham.col, ham.row, ham.data]))
            ds.attrs['n_fock'] = len(states)

# ----
# Loop through spectra.
spectrum_file = (
    param['working_directory'] +
    '/multi_spectrum_' +
    param['gauge_particles'] + '_' +
    ('' if ws is None else ws+'_') +
    size_tag(param['L']) +
    '.dat'
)

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

    # Output.
    with hdf.File(spectrum_file, 'a' if i else 'w') as f:
        f.attrs['lambdas'] = lambdas
        ds = f.create_dataset('lam_{:6f}'.format(l), data=spectra[l])
        ds.attrs['lambda'] = l


if args.notify:
    from gauss_lattice.lr_notify import push_message
    import os, subprocess
    host = subprocess.check_output(['hostname']).strip().decode('UTF-8'),
    push_message(f'Multi-parameter run  @{host[0]} is done!')
