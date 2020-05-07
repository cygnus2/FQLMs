""" ----------------------------------------------------------------------------

    hamiltonian_builder.py - LR, May 2020

    Constructs a hamiltonian for the gauss lattice.

---------------------------------------------------------------------------- """
import numpy as np
from bisect import bisect_left
from .hamiltonian import Hamiltonian
from .bit_magic import set_bits, sum_occupancies
from copy import copy

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

        # Set the type of partilces we're considering (only changes the sign
        # factors in the Hamiltonian).
        gp = param.get('gauge_particles')
        self.fermions = gp == 'fermions'

        # Pre-compute all plaquette indicies to save some time.
        self.plaquettes = self.get_plaquette_list()


        # Prepare for operator application.
        self.op_u = [False, False, True, True]
        self.op_u_dagger = [True, True, False, False]


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
            plaquettes.append(ind + [set_bits(ind)])

            # In 3D, we have two additional plaquettes.
            if self.d == 3:
                # yz plane
                j = self.shift_index(self.shift_index(n, 2), 1) # shifted by Sy and Sz
                vn_yz = self.get_vertex_links(j)
                ind=[vn[2], vn_yz[5], vn_yz[3], vn[4]]
                plaquettes.append(ind + [set_bits(ind)])

                # xz plane.
                j = self.shift_index(self.shift_index(n, 0), 2) # shifted by Sx and Sz
                vn_xz = self.get_vertex_links(j)
                ind = [vn[0], vn_xz[5], vn_xz[1], vn[4]]
                plaquettes.append(ind + [set_bits(ind)])

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

            # First apply the U term.
            # new_state, sign = self.apply_u_dagger(state, p)
            new_state, sign = self.apply_plaquette_operator(state, p, mask=self.op_u)

            # If U term was not successful, try the U^dagger term.
            # (the order could have been switched - there's always only one
            # possibility for overlap to be generated)
            if not new_state:
                # new_state, sign = self.apply_u(state, p)
                new_state, sign = self.apply_plaquette_operator(state, p, mask=self.op_u_dagger)

            if new_state:
                states.append([
                    n_state,
                    self.state_to_index(new_state),
                    sign if self.fermions else 1
                ])
        return states


    def apply_plaquette_operator(self, state, p, mask):
        """ Applies the U operator

                c1+ c2+ c3 c4

            or the U^dagger term

                c1 c2 c3+ c4+

            to a given plaquette in a given state (link configuration starting
            with x-mu link and going counter-clockwise).
        """
        a, b, n = max(p[:-1]), min(p[:-1]), 0
        new_state = copy(state)
        for k in range(4):
            m = 1 << p[k]
            if bool(new_state & m) == mask[k]:
                if self.fermions:
                    n += sum_occupancies(a, p[k], new_state)
                new_state = copy(new_state^m)
            else:
                return 0, 0
        return new_state, (-1)**n


    def apply_u(self, state, p):
        """ Applies the U operator term

                c_1^+ c_2^+ c_3 c_4

            to a given plaquette in a given state (link configuration starting
            with x-mu link and going counter-clockwise).

            Logic: If links 1&2 are occupied and links 2&3 are free, then we can
                   apply U operation.
        """
        a, b = max(p[:-1]), min(p[:-1])

        m = 1 << p[0]
        if not state & m:
            n = sum_occupancies(a, p[0], state)
            new_state = copy(state^m)

            m = 1 << p[1]
            if not new_state & m:
                n += sum_occupancies(a, p[1], new_state)
                new_state = copy(new_state^m)

                m = 1 << p[2]
                if new_state & m:
                    n += sum_occupancies(a, p[2], new_state)
                    new_state = copy(new_state^m)

                    m = 1 << p[3]
                    if new_state & m:
                        n += sum_occupancies(a, p[3], new_state)
                        new_state = copy(new_state^m)
                        return new_state, (-1)**n
        return 0, 0


    def apply_u_dagger(self, state, p):
        """ Applies the U^dagger opator term
            to a given plaquete in a given state.
            Logic: If links 1&2 are free and links 2&3 are occupied, then we can
                   apply the U^dagger operation.
        """
        a, b = max(p[:-1]), min(p[:-1])

        m = 1 << p[0]
        if state & m:
            n = sum_occupancies(a, p[0], state)
            new_state = copy(state^m)

            m = 1 << p[1]
            if new_state & m:
                n += sum_occupancies(a, p[1], new_state)
                new_state = copy(new_state^m)

                m = 1 << p[2]
                if not new_state & m:
                    n += sum_occupancies(a, p[2], new_state)
                    new_state = copy(new_state^m)

                    m = 1 << p[3]
                    if not new_state & m:
                        n += sum_occupancies(a, p[3], new_state)
                        new_state = copy(new_state^m)
                        return new_state, (-1)**n
        return 0, 0



    def index_to_state(self, n):
        """ Maps the n-th state in the Fock basis to it's bit string.
        """
        return self.lookup_table[n]


    def state_to_index(self, state):
        """ Maps the bit string to the basis index. Lookup is based on bisection
            which requires an ordered list (ascending, I think, is mandatory).
        """
        return bisect_left(self.lookup_table, state)
