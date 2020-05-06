""" ----------------------------------------------------------------------------

    hamiltonian_builder.py - LR, May 2020

    Constructs a hamiltonian for the gauss lattice.

---------------------------------------------------------------------------- """
import numpy as np
import datetime as dt
from bisect import bisect_left
from hamiltonian import Hamiltonian


class HamiltonianBuilder(object):
    """ Constructs the Hamiltonian in a general single-particle basis.
    """

    def __init__(self, param, states):
        """
        """
        self.L = param['L']

        # Shift operators to get the link indices.
        self.S = [1]
        for l in self.L:
            self.S.append(l*self.S[-1])


        # Set up the lookup table.
        self.generate_lookup_table()
        self.n_fock = len(self.lookup_table)

        # Pre-compute all plaquette indicies to save some time.
        self.plaquettes = self.get_plaquette_list()


    def generate_lookup_table(self):
        """ Produces the ordered lookup table for the states.
        """
        # raise NotImplementedError("sorry")
        self.lookup_table = [1]


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
        # raise NotImplementedError("sorry")
        return []


    def construct(self, timing=True):
        """ Actually builds the Hamiltonian and returns a Hamiltonian object
            ready to be diagonalized.
        """
        if timing:
            print('Starting construction of Hamiltonian matrix.')
            start = dt.datetime.now()

        # Loop through all Fock states and create the overlap matrix. First step:
        # do it naively (with some doubled work). Then try to improve on that (by
        # using, e.g., Hermiticity).
        all_entries = []
        for n in range(self.n_fock):
            all_entries += self.do_single_state(n)

        #DEV: take out!
        all_entries = [[0, 0, 1]]

        # Make a sparse matrix out ot this -although pretty plain, this can handle
        # reasonably sized lists of indices (will do fo now).
        row, col, data = zip(*all_entries)
        ham = Hamiltonian(data, row, col, shape=(self.n_fock, self.n_fock))

        if timing:
            print('Done with construction of Hamiltonian matrix.')
            duration = dt.datetime.now() - start
            print(f'[runtime was {duration}]')
        return ham


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
            pass
            # U operator.

            # U^\dagger operator.

        return states


    def index_to_state(self, n):
        """ Maps the n-th state in the Fock basis to it's bit string.
        """
        return self.lookup_table[n]

    def state_to_index(self, state):
        """ Maps the bit string to the basis index. Lookup is based on bisection.
        """
        return bisect_left(self.lookup_table, state)