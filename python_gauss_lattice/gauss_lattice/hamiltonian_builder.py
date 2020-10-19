""" ----------------------------------------------------------------------------

    hamiltonian_builder.py - LR, May 2020

    Constructs a hamiltonian for the gauss lattice.

---------------------------------------------------------------------------- """
import numpy as np
from bisect import bisect_left
from .hamiltonian import GaussLatticeHamiltonian
from .hamiltonian_builder_methods import do_single_state
from .bit_magic import set_bits, sum_occupancies
from .aux import timestamp
from copy import copy
from tqdm import tqdm as tbar
from multiprocessing import Pool
from itertools import product


class HamiltonianBuilder(object):
    """ Constructs the Hamiltonian in a general single-particle basis.
    """

    def __init__(self, param, states, logger=None, silent=False):
        self.L = param['L']
        self.d = len(param['L'])

        # Shift operators to get the link indices.
        self.S = [1]
        for l in self.L:
            self.S.append(l*self.S[-1])

        # Flag for big int storage.
        self.nb = self.S[-1]*self.d
        self.big_int = self.nb > 60

        # Some I/O business.
        self.silent = silent
        self.logger = logger

        # Set up the lookup table, which is merely ordering the states such that
        # the inverse lookup can be done efficiently with bisection.
        self.lookup_table = sorted(states)
        self.n_fock = len(self.lookup_table)
        if not self.silent:
            self._log(f'Setting up the Hamiltonian with {self.n_fock} Fock states.')

        # Pre-compute all plaquette indicies to save some time.
        self.plaquettes = self.get_plaquette_list()



    def _log(self, msg):
        """ To produce a logfile (mainly debugging purposes).
        """
        if self.logger:
            self.logger.info(timestamp() + ' ' + msg)
        else:
            print(timestamp() + ' ' + msg)


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

            (ASCII art shamelessly stolen from Debasish)

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


    def construct(self, n_threads=1, progress_bar=False):
        """ Actually builds the Hamiltonian and returns a Hamiltonian object
            ready to be diagonalized.
        """
        # Loop through all Fock states and create the overlap matrix. First step:
        # do it naively (with some doubled work). Then try to improve on that (by
        # using, e.g., Hermiticity).

        self._log(f'Working with {n_threads} threads.')
        all_entries = []
        if n_threads == 1:
            if progress_bar:
                for s in tbar(self.lookup_table):
                    all_entries += [do_single_state((s, self.plaquettes))]
            else:
                for s in self.lookup_table:
                    all_entries += [do_single_state((s, self.plaquettes))]

        else:
            with Pool(n_threads) as pool:
                all_entries = pool.map(do_single_state, product(self.lookup_table, [self.plaquettes]))

        # Make a sparse matrix out ot this -although pretty plain, this can handle
        # reasonably sized lists of indices (will do fo now).
        icol, irow, idata = [], [], []
        if len(all_entries):
            for line in all_entries:
                if len(line):
                    row, col, data = zip(*line)

                    # Convert the columns to the proper format.
                    for k in range(len(row)):
                        c = self.state_to_index(col[k])
                        if c < self.n_fock:
                            icol.append(c)
                            irow.append(self.state_to_index(row[k]))
                            idata.append(data[k])

        if not self.silent:
            self._log("# of nonzero entries: " + str(len(idata)))
        return GaussLatticeHamiltonian(idata, irow, icol, n_fock=self.n_fock)


    def index_to_state(self, n):
        """ Maps the n-th state in the Fock basis to it's bit string.

            Note:
                - likely deprecated, here for historic purposes.
        """
        return self.lookup_table[n]


    def state_to_index(self, state):
        """ Maps the bit string to the basis index. Lookup is based on bisection
            which requires an ordered list (ascending, I think, is mandatory).
        """
        return bisect_left(self.lookup_table, state)
