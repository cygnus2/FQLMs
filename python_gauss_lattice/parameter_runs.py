""" ----------------------------------------------------------------------------

    parameter_runs.py - LR, May 2020

    Scans through parameter space starting from a constructed Hamiltonian.

---------------------------------------------------------------------------- """
from gauss_lattice import GaussLatticeHamiltonian, HamiltonianBuilder, GaussLattice
from gauss_lattice.aux import size_tag, timeit, read_all_states, file_tag
import numpy as np
import h5py as hdf

@timeit
def hamiltonian_diagonalization(ham, **kwargs):
    """ This method just exists to be able to time the routine efficiently.
    """
    return ham.diagonalize(**kwargs)

# ---

# This sets the parameters fof the calculation (everything else is fixed).
param = {
    'L' : [4,4],
    'J' : -1.0,
    'gauge_particles' : 'bosons',

    'ev_type' : 'BE',
    'n_eigenvalues' : 100,
    'full_diag' : False
}
lambdas = np.linspace(-0.1, 0.1, 100)

# ---
# Retrieves the Hamiltonian and if necessary, constructs it.

input_file = 'output/hamiltonian_' + size_tag(param['L']) + '.npz'
try:
    ham = GaussLatticeHamiltonian.from_scipy_dump(input_file)
except FileNotFoundError:
    try:
        states = read_all_states(param['L'])
    except FileNotFoundError:
        glatt = GaussLattice(param['L'], state_file=file_tag(L, filetype='hdf5'))
        glatt.find_states()
        states = read_all_states(param['L'])
    builder = HamiltonianBuilder(param, states)
    ham = builder.construct()
    ham.store_hamiltonian(input_file)

# ----

spectra = {}
for l in lambdas:
    # Diagonalization.
    spectra[l] = hamiltonian_diagonalization(ham,
        full_diag = param.get('full_diag'),
        J = param['J'],
        lam = l,
        gauge_particles = param['gauge_particles'],
        n_eigenvalues = max(1, min(param['n_eigenvalues'], ham.n_fock//2)),
        which = param['ev_type']
    )

with hdf.File('output/multi_spectrum_'+size_tag(param['L']) + '.hdf5', 'w') as f:
    f.attrs['lambdas'] = lambdas
    for l in lambdas:
        f.create_dataset('lam_{:6f}'.format(l), data=spectra[l])
