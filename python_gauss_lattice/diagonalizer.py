""" ----------------------------------------------------------------------------

    diagonalizer.py - LR, May 2020

    Run script for diagonalizing the Hamiltonian for the gauss lattice.

---------------------------------------------------------------------------- """
from gauss_lattice import HamiltonianBuilder
from gauss_lattice.aux import size_tag, timeit, read_winding_sector, read_all_states
import numpy as np

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

# This sets the parameters fof the calculation (everything else is fixed).
param = {
    'L' : [2,2,2]
    # 'winding_sector' : (0,0),
}
if param.get('winding_sector'):
    states = read_winding_sector(param['L'], param['winding_sector'])
else:
    states = read_all_states(param['L'])

# Set up the builder object & construct the Hamiltonian.
builder = HamiltonianBuilder(param, states=states)
ham = hamiltonian_construction(builder)

# Just to have it.
ham.store_hamiltonian('output/hamiltonian_'+size_tag(param['L'])+'.npz')

# ----------------------
# Diagonalization.

particles = 'bosons'
full_diag = False
n_eigenvalues = min(100, builder.n_fock//2)

eigenvalues = hamiltonian_diagonalization(
    ham,
    gauge_particles='bosons',
    lam=1,
    full_diag=full_diag,
    n_eigenvalues=n_eigenvalues,
    which='BE'
)

# Some I/O.
filename = 'spectrum_' + particles + '_' + size_tag(param['L']) + '.dat'
if full_diag:
    filename = 'FULL_' + filename
ham.store_results(filename='output/'+filename)
