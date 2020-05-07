""" ----------------------------------------------------------------------------

    hamiltonian_builder.py - LR, May 2020

    Constructs a hamiltonian for the gauss lattice.

---------------------------------------------------------------------------- """
import numpy as np
from bisect import bisect_left
from .hamiltonian import Hamiltonian


class HamiltonianBuilder(object):
    """ Constructs the Hamiltonian in a general single-particle basis.
    """

    def __init__(self, param, states):
        """
        """
        self.L = param['L']
        self.d = len(param['L'])

        # Shift operators to get the link indices.
        self.S = [1]
        for l in self.L:
            self.S.append(l*self.S[-1])

        # Set up the lookup table, which is merely ordering the states such that
        # the inverse lookup can be done efficiently with bisection.
        self.lookup_table = sorted(states)
        self.n_fock = len(self.lookup_table)
        print(f'Setting up the Hamiltonian with {self.n_fock} Fock states.')

        # Pre-compute all plaquette indicies to save some time.
        self.plaquettes = self.get_plaquette_list()


    def shift_index(self, i, d):
        """ Shifts an grid index into direction d under consideration of  PBC.
            Works only for positive shifts (in positive coordinate axes) but
            this is also the only use case.
        """
        L, S = self.L, self.S
        d += 1

        n = i + S[d-1]
        if n % S[d] < S[d-1]:
            n -= S[d] #- S[d-1]
        return n


    def get_plaquette_list(self):
        """ Produces a list of plaquettes that represent the lattice. This is a
            list of plaquettes represented as

                plaquette = [link_1, link_2, link_3, link_4]

            where link link_1 and link_2 as well as link_3 and link_4 are
            connected through the hopping like

                                      link3
                                    O--------O
                                    |        |
                              link4 |        | link2
                                    |        |
                                    O--------O
                                      link1

            (ASCII art stolen from Debasish)

            This list will later be used when we loop through the states to
            construct the Hamiltonian.
        """
        S = self.S

        # Generates a bit mask for later use when we flip plaquettes.
        make_mask = lambda ilist: sum([1<<k for k in ilist])

        # Find plaquettes by looping over all grid points.
        plaquettes = []
        for n in range(S[-1]):

            # Contents of a vertex:
            # 0  1  2  3  4  5 ...
            # x -x  y -y  z -z ...
            vn = self.get_vertex_links(n) # on site

            # xy plane.
            j = self.shift_index(self.shift_index(n, 0), 1) # shifted by Sx and Sy
            vn_xy = self.get_vertex_links(j)
            ind = [vn[0], vn_xy[3], vn_xy[1], vn[2]]
            plaquettes.append(ind + [make_mask(ind)])

            # In 3D, we have two additional plaquettes.
            if self.d == 3:
                # yz plane
                j = self.shift_index(self.shift_index(n, 2), 1) # shifted by Sy and Sz
                vn_yz = self.get_vertex_links(j)
                ind=[vn[2], vn_yz[5], vn_yz[3], vn[4]]
                plaquettes.append(ind + [make_mask(ind)])

                # xz plane.
                j = self.shift_index(self.shift_index(n, 0), 2) # shifted by Sx and Sz
                vn_xz = self.get_vertex_links(j)
                ind = [vn[0], vn_xz[5], vn_xz[1], vn[4]]
                plaquettes.append(ind + [make_mask(ind)])

        # Check if the right amount of plaquettes was found and if so, return
        # the list.
        assert len(plaquettes) == (2**(self.d-1) -1) * S[-1]
        return plaquettes



    def get_vertex_links(self, i):
        """ Returns the indices [+x, -x, +y, -y, ...] of the links for the i-th
            vertex *on the full lattice* (not the sublattice) under consideration
            of periodic boundary conditions.

            Note: works only for L^d lattices for now.
        """
        # Shorthand to avoid self all the time.
        L, d, S = self.L, self.d, self.S

        # Index in bit string.
        j = d*i

        # Loop through the dimensions and add the indicies of + and - directions.
        ind = []
        for k in range(1,d+1):
            # Step forward is always 'on site'.
            ind.append(j+k-1)

            # Sep backward in k direction.
            l = j + k -1 - d*S[k-1]
            if i%S[k] < S[k-1]:
                l = l + d*S[k]
            ind.append(l)
        return ind


    def construct(self):
        """ Actually builds the Hamiltonian and returns a Hamiltonian object
            ready to be diagonalized.
        """
        # Loop through all Fock states and create the overlap matrix. First step:
        # do it naively (with some doubled work). Then try to improve on that (by
        # using, e.g., Hermiticity).
        all_entries = []
        for n in range(self.n_fock):
            all_entries += self.do_single_state(n)

        # Make a sparse matrix out ot this -although pretty plain, this can handle
        # reasonably sized lists of indices (will do fo now).
        row, col, data = zip(*all_entries)
        return Hamiltonian(data, row, col, shape=(self.n_fock, self.n_fock))


    def do_single_state(self, n_state):
        """ Constructs matrix-elements for a single Fock state.
        """
        # First get bit representation.
        state = self.index_to_state(n_state)

        # This loop applies all the plaquette terms of the Hamiltonian. The
        # strategy is as follows:
        #   1 - loop through all the plaquettes
        #   2 - for every plaquette, apply both operators (U and U^dagger) only
        #       if the plaquette is "flipable" - this condition will be checked
        #       on the fly.
        #   3 - inverse lookup (based on bisection) to determine the index of
        #       the newly generated states.
        #   4 - return the list to append to the sparse matrix representation.
        states = []
        for p in self.plaquettes:
            flip = False

            # If links 1&2 are occupied and links 2&3 are free, then we can apply
            # the U operation.
            if state & (1<<p[0]):
                if state & (1<< p[1]):
                    if not state & (1 << p[2]):
                        if not state & (1 << p[3]):
                            # Here we flip the spins that are in the mask.
                            flip = True

            # If links 1&2 are free and links 2&3 are occupied, then we can apply
            # the U^dagger operation.
            if not state & (1<<p[0]):
                if not state &( 1<< p[1]):
                    if state & (1 << p[2]):
                        if state & (1 << p[3]):
                            # Here we flip the spins that are in the mask.
                            flip = True

            if flip:
                new_state = state^p[-1]

                sign = 1
                states.append([
                    n_state,
                    self.state_to_index(new_state),
                    sign
                ])

        return states


    def index_to_state(self, n):
        """ Maps the n-th state in the Fock basis to it's bit string.
        """
        return self.lookup_table[n]


    def state_to_index(self, state):
        """ Maps the bit string to the basis index. Lookup is based on bisection
            which requires an ordered list (ascending, I think, is mandatory).
        """
        return bisect_left(self.lookup_table, state)
