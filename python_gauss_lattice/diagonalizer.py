""" ----------------------------------------------------------------------------

    diagonalizer.py - LR, May 2020

    Run script for diagonalizing the Hamiltonian for the gauss lattice.

---------------------------------------------------------------------------- """
from gauss_lattice import HamiltonianBuilder
from gauss_lattice.aux import size_tag, timeit, read_winding_sector, read_all_states, load_config
import numpy as np
import argparse


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

# ------------------------------------------------------------------------------
# Hamiltonian setup.
if param.get('winding_sector'):
    states = read_winding_sector(param['L'], param['winding_sector'], basedir=param['working_directory'])
else:
    states = read_all_states(param['L'], basedir=param['working_directory'])

# Set up the builder object & construct the Hamiltonian.
builder = HamiltonianBuilder(param, states=states)
ham = hamiltonian_construction(builder)

# Just to have it.
ham.store_hamiltonian(param['working_directory'] + 'hamiltonian_'+size_tag(param['L'])+'.npz')


# ------------------------------------------------------------------------------
# Diagonalization.
results = hamiltonian_diagonalization(ham,
    gauge_particles = param['gauge_particles'],
    J = param['J'],
    lam = param['lambda'],
    full_diag = param.get('full_diag'),
    n_eigenvalues = min(param['n_eigenvalues'], builder.n_fock//2),
    which = param['ev_type'],
    compute_eigenstates=bool(param.get('compute_eigenstates'))
)
if param['compute_eigenstates']:
    eigenvalues, eigenstates = results
else:
    eivenvalues = results

# ------------------------------------------------------------------------------
# Some I/O.
filename = 'spectrum_' + param['gauge_particles'] + '_' + size_tag(param['L']) + '_lam{:.2f}'.format(param['lambda']) + '.dat'
if param.get('full_diag'):
    filename = 'FULL_' + filename
ham.store_results(filename=param['working_directory']+filename, store_eigenvalues=True)
