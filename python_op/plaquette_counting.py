import numpy as np
import h5py as hdf

import sys
sys.path.append('../python_gauss_lattice/')
from gauss_lattice.hamiltonian_builder_methods import apply_u_dagger, apply_u
from gauss_lattice.aux_stuff import timestamp, read_winding_sector
from gauss_lattice import HamiltonianBuilder


def count_flippable_plaquettes(state, plaquettes):
    """ Finds the flippable plaquettes and returns them.
    """
    c_1, c_2 = 0, 0
    for p in plaquettes:
        if apply_u_dagger(state, p, sign=False)[0]:
            c_1 += 1
        elif apply_u(state, p, sign=False)[0]:
            c_2 += 1
    return c_1, c_2


# Read the GL states for a given winding sector.
L = [2,2,4]
winding_states, _ = read_winding_sector(L, [0,0,0], basedir='../python_data/local_state_storage/')
winding_states = np.sort(winding_states)

# Count.
builder = HamiltonianBuilder({'L':L}, [], silent=True)
n_states = len(winding_states)
plaquettes = builder.get_plaquette_list(separate_lists=True)

# Compute observables.
obs = np.zeros(shape=(n_states,6))
for k in range(n_states):
    for j, p in enumerate(plaquettes):
        obs[k,2*j:2+2*j] = count_flippable_plaquettes(winding_states[k], p)

    # Show progress.
    if not (k+1)%(1e4):
        print(f'{timestamp()} {k+1} / {n_states} states done')

np.save('counted_plaquettes', obs)
