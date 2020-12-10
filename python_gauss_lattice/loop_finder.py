""" ----------------------------------------------------------------------------

    loop_finder.py - LR, May 2020

    Tries to find negatively contributing worldlines. This is quite memory
    inefficient and should only be done for small systems (in 3D, only 2x2x2,
    I guess).

---------------------------------------------------------------------------- """
import numpy as np
from scipy.sparse import load_npz
from gauss_lattice.aux_stuff import size_tag, read_all_states
from gauss_lattice import HamiltonianBuilder
from gauss_lattice.bit_magic import count_particles

def print_path(path):
    print(' >>> '.join(map(str, path)))

param = {
    'L' : [2,2,2],
}
input_file = 'output/hamiltonian_' + size_tag(param['L']) + '.npz'
ham = load_npz(input_file)

# Read the GSL states and produce a builder object, this is convenient for later.
gsl_states = read_all_states(param['L'])
builder = HamiltonianBuilder(param, states=gsl_states)
blen = builder.S[-1]*3

# Here we need to convert the sparse matrix into a list of transitions, since
# this is easier to use for the recursion.
n_fock = ham.shape[0]
transitions = {k : [] for k in range(n_fock)}
for i, r in enumerate(ham.row):
    transitions[r].append((ham.col[i], int(ham.data[i])))

def loop_state(path, ml):
    if len(path) >= ml:
        return []
    for i, n in enumerate(path[:-1]):
        if abs(path[-1]) == abs(n):
            if len(path)-i > 3:
                if np.sign(path[i]) != np.sign(path[-1]):
                    return [path[i:]]
            return []
    res = []
    for next in transitions[abs(path[-1])]:
        res += loop_state(path+[np.sign(path[-1])*next[1]*next[0]], ml)
    return res


# Actually finding the loops (without checking anything).
n_state = 4815
loops = loop_state([n_state], 8)[:1]

# n_particles = 12
# for state in gsl_states:
#     if count_particles(state) == n_particles:
#         n_state = builder.state_to_index(state)
#         print('{:024b}'.format(state))
#         loops = loop_state([n_state], 8)[:1]
#         if len(loops):
#             break

# Here we prepare the loops for efficient later use.
#   1) which indicies changed?
#   2) which plaquette
#   3) and what operator was that? (dagger or not)
#
# This is a small abuse of Code, but it's convenient to use the HAmiltonianBuilder
# for that.
bit_loops = [list(map(lambda x: builder.index_to_state(abs(x)), loop)) for loop in loops]
ind_loops = [list(map(lambda x: [k for k in range(blen) if (x>>k)&1], bl)) for bl in bit_loops]

def find_plaquette(ilist, plaquettes):
    for p in plaquettes:
        if set(ilist) == set(p[:-1]):
            return p[:-1]

# Get the list of affected plaquettes & find the type of operator that was applied..
plaquette_strings, operator_strings = [], []
for bl in bit_loops:
    plaquettes, operators = [], []
    for i in range(len(bl)-1):
        diff = bl[i] ^ bl[i+1]
        ilist = [k for k in range(blen) if (diff>>k)&1]
        plaquette = find_plaquette(ilist, builder.plaquettes)
        plaquettes.append(plaquette)
        if (bl[i] >> plaquette[0])&1 :
            operators.append('U+')
        else:
            operators.append('U ')
    operator_strings.append(operators)
    plaquette_strings.append(plaquettes)

for k, loop in enumerate(loops):
    # path = ''
    # for j in range(len(loop)-1):
    #     path += '{:s}[{:2d},{:2d},{:2d},{:2d}] --> '.format(operator_strings[k][j], *plaquette_strings[k][j])
    # print(path[:-5])

    # More explicit.
    for j in range(len(loop)-1):
        st = '[ {:d} --> {:d} ] '.format(loop[j], loop[j+1])
        op = '{:s}[{:2d},{:2d},{:2d},{:2d}] '.format(operator_strings[k][j], *plaquette_strings[k][j])
        state = '| '+ ('{:2d}, '*len(ind_loops[k][j])).format(*ind_loops[k][j][::-1])[:-2] + '  >'
        print(st + op + ' ' + state)
