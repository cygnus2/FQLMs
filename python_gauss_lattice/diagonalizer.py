""" ----------------------------------------------------------------------------

    diagonalizer.py - LR, May 2020

    Run script for diagonalizing the Hamiltonian for the gauss lattice.

---------------------------------------------------------------------------- """
from gauss_lattice import HamiltonianBuilder, GaussLatticeHamiltonian
from gauss_lattice.aux import size_tag, timeit, read_winding_sector, read_all_states, load_config
import numpy as np
import argparse, logging
import h5py as hdf


@timeit(logger=None)
def hamiltonian_construction(builder):
    """ This method just exists to be able to time the routine efficiently.
    """
    return builder.construct()

@timeit(logger=None)
def hamiltonian_diagonalization(ham, **kwargs):
    """ This method just exists to be able to time the routine efficiently.
    """
    return ham.diagonalize(**kwargs)

# ------------------------------------------------------------------------------
# Input handling.
parser = argparse.ArgumentParser(description="Python gauss lattice diagonalizer (single lambda).")
parser.add_argument('-i', metavar='', type=str, default=None, help='YAML style input file.')
args = parser.parse_args()

# Get parameters.
param = load_config(args.i)
hamiltonian_file = param['working_directory'] + '/SEQUENTIAL_hamiltonian_' + size_tag(param['L']) + '.hdf5'

# Set up a logger with a handler for the terminal output.
logger = logging.getLogger('diagonalization logger')
logger.addHandler(logging.StreamHandler())
logger.addHandler(logging.FileHandler(param['working_directory'] + "/" + param['logfile'], mode='w'))
logger.setLevel(logging.DEBUG)

# ------------------------------------------------------------------------------
# Hamiltonian setup.
if param.get('winding_sector'):
    states, ws = read_winding_sector(param['L'], param['winding_sector'], basedir=param['working_directory'])
else:
    states = read_all_states(param['L'], basedir=param['working_directory'])
    ws = None

# Either get an existing Hamiltonian matrix or produce it from scratch.
# Note: only specific winding sectors can be rwead, the full hamiltonian cannot be
# read from file. This is OK, since the real benefit only exists at large systems
# (2x2x4 mainly), whereas small systems can easily be constructed with little effort.
try:
    with hdf.File(hamiltonian_file, 'r') as f:
        mat = f[ws][...]
    if mat.shape[-1]:
        ham = GaussLatticeHamiltonian(mat[2,:], mat[1,:], mat[0,:], len(states))
    else:
        ham = GaussLatticeHamiltonian([], [], [], 1)
    logger.info('Read Hamiltonian from provided file.')

    # ham = GaussLatticeHamiltonian(data, row, col, nfock)
except (KeyError, OSError):
    # Set up the builder object & construct the Hamiltonian.
    builder = HamiltonianBuilder(param, states=states, logger=logger)
    ham = hamiltonian_construction(builder)

    # If specified, store the Hamiltonian for later use.
    if param['store_hamiltonian']:
        with hdf.File(hamiltonian_file, 'a' if i else 'w') as f:
            ds = f.create_dataset(ws, data= np.array([ham.col, ham.row, ham.data]))
            ds.attrs['n_fock'] = len(states)


# ------------------------------------------------------------------------------
# Diagonalization.
n_ev = min(param['n_eigenvalues'], ham.n_fock//2)
logger.info('Starting to diagonalize, computing {:d} eigenvalues.'.format(n_ev))
results = hamiltonian_diagonalization(ham,
    gauge_particles = param['gauge_particles'],
    J = param['J'],
    lam = param['lambda'],
    full_diag = param.get('full_diag'),
    n_eigenvalues = n_ev,
    which = param['ev_type'],
    compute_eigenstates=bool(param.get('compute_eigenstates'))
)
if param['compute_eigenstates']:
    eigenvalues, eigenstates = results
else:
    eivenvalues = results

# ------------------------------------------------------------------------------
# Some I/O.
filename = (
    'spectrum_' +
    param['gauge_particles'] + '_' +
    ('' if ws is None else ws+'_') +
    size_tag(param['L']) +
    '_lam{:.2f}'.format(param['lambda']) + '.dat'
)
if param.get('full_diag'):
    filename = 'FULL_' + filename
ham.store_results(filename=param['working_directory']+filename, store_eigenvalues=True)
