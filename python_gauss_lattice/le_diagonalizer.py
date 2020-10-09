""" ----------------------------------------------------------------------------

    diagonalizer.py - LR, May 2020

    Run script for diagonalizing the Hamiltonian for the gauss lattice.

---------------------------------------------------------------------------- """
from gauss_lattice import LowEnergyHamiltonianBuilder
from gauss_lattice.aux import size_tag, timeit, load_config
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

# Set up a logger with a handler for the terminal output.
logger = logging.getLogger('diagonalization logger')
logger.addHandler(logging.StreamHandler())
logger.addHandler(logging.FileHandler(param['working_directory'] + "/" + param['logfile'], mode='w'))
logger.setLevel(logging.DEBUG)

# ------------------------------------------------------------------------------
# Hamiltonian setup.
builder = LowEnergyHamiltonianBuilder(param, logger=logger)
base_lattices = {
    (2,2,2) : [3816540, 3872106, 5421780, 5678001, 7542990, 7743645,
                9033570, 9234225, 11099214, 11355435, 12905109, 12960675],
    (2,2,4) : [64030919769180, 64963162609002, 90962379586260, 95261054903217,
                126550380058830, 129916812535965, 151558164174690, 154924596651825,
                186213921807438, 190512597124395, 216511814101653, 217444056941475]
}
builder.find_le_states(base_lattices[tuple(param['L'])], param['excitation_levels'])
ham = hamiltonian_construction(builder)

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
    'le_spectrum_' +
    param['gauge_particles'] + '_' +
    "ex{:d}_".format(param["excitation_levels"]) +
    size_tag(param['L']) +
    '_lam{:.2f}'.format(param['lambda']) + '.dat'
)
if param.get('full_diag'):
    filename = 'FULL_' + filename
ham.store_results(filename=param['working_directory']+filename, store_eigenvalues=True)
