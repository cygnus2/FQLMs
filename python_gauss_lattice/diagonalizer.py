""" ----------------------------------------------------------------------------

    diagonalizer.py - LR, May 2020

    Run script for diagonalizing the Hamiltonian for the gauss lattice.

---------------------------------------------------------------------------- """
from gauss_lattice import HamiltonianBuilder
from gauss_lattice.aux import size_tag, timeit, read_winding_sector, read_all_states
import numpy as np
import yaml, argparse


@timeit
def hamiltonian_construction(builder):
    """ This method just exists to be able to time the routine efficiently.
    """
    return builder.construct()

@timeit
def hamiltonian_diagonalization(ham, **kwargs):
    """ This method just exists to be able to time the routine efficiently.
    """
    return ham.diagonalize(**kwargs)

# ------------------------------------------------------------------------------
# Input handling.

parser = argparse.ArgumentParser(description="Python gauss lattice diagonalizer.")
parser.add_argument('-i', metavar='', type=str, default=None, help='YAML style input file.')
args = parser.parse_args()

if args.i is not None:
    # Loads parameters from file.
    with open(args.i, 'r') as f:
        try:
            param = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            print(exc)
else:
    sys.exit('fatal: no input file specified')


# ------------------------------------------------------------------------------
# Hamiltonian setup.
if param.get('winding_sector'):
    states = read_winding_sector(param['L'], param['winding_sector'])
else:
    states = read_all_states(param['L'])

# Set up the builder object & construct the Hamiltonian.
builder = HamiltonianBuilder(param, states=states)
ham = hamiltonian_construction(builder)

# Just to have it.
ham.store_hamiltonian('output/hamiltonian_'+size_tag(param['L'])+'.npz')


# ------------------------------------------------------------------------------
# Diagonalization.
eigenvalues = hamiltonian_diagonalization(ham,
    gauge_particles = param['gauge_particles'],
    J = param['J'],
    lam = param['lambda'],
    full_diag = param['full_diag'],
    n_eigenvalues = min(param['n_eivenvalues'], builder.n_fock//2),
    which = param['ev_type']
)


# ------------------------------------------------------------------------------
# Some I/O.
filename = 'spectrum_' + param['gauge_particles'] + '_' + size_tag(param['L']) + '_lam{:.2f}'.format(param['lambda']) + '.dat'
if param['full_diag']:
    filename = 'FULL_' + filename
ham.store_results(filename='output/'+filename)
