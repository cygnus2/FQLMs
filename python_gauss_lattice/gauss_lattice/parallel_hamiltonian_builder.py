""" ----------------------------------------------------------------------------

    parallel_hamiltonian_builder.py - LR, October 2020

    More efficient parallel version of the Hamiltonian builder.

---------------------------------------------------------------------------- """
from gauss_lattice.hamiltonian_builder_methods import apply_u, apply_u_dagger
from gauss_lattice.hamiltonian_builder import HamiltonianBuilder
from gauss_lattice import GaussLatticeHamiltonian
from multiprocessing import Pool

class ParallelHamiltonianBuilder(HamiltonianBuilder):
    """ Constructs the Hamiltonian in a general single-particle basis.
    """
    plaquettes = None
    inv_lookuptable = None

    def __init__(self, *args, **kwargs):
        HamiltonianBuilder.__init__(self, *args, **kwargs)
        ParallelHamiltonianBuilder.set_plaquettes(self.plaquettes)
        ParallelHamiltonianBuilder.set_inv_lookuptable(self.lookup_table)

    @staticmethod
    def set_plaquettes(plaquettes):
        ParallelHamiltonianBuilder.plaquettes = plaquettes

    @staticmethod
    def set_inv_lookuptable(lookup_table):
        ParallelHamiltonianBuilder.inv_lookuptable = {v:k for k,v in enumerate(lookup_table)}

    @staticmethod
    def do_single_state(state, sign=True):
        """ Tries to flip all plaquettes in a single state.
        """
        states = []
        for p in ParallelHamiltonianBuilder.plaquettes:
            # First apply the U term.
            new_state, s = apply_u_dagger(state, p, sign=sign)

            # If U term was not successful, try the U^dagger term.
            # (the order could have been switched - there's always only one
            # possibility for overlap to be generated)
            if not new_state:
                new_state, s = apply_u(state, p, sign=sign)

            if new_state:
                c = ParallelHamiltonianBuilder.inv_lookuptable.get(new_state)
                if c:
                    states.append([state, c, s])

        return states


    def construct(self, n_threads=1):
        """ Actually builds the Hamiltonian and returns a Hamiltonian object
            ready to be diagonalized.
        """
        # Loop through all Fock states and create the overlap matrix. First step:
        # do it naively (with some doubled work). Then try to improve on that (by
        # using, e.g., Hermiticity).

        self._log(f'Working with {n_threads} threads.')
        all_entries = []
        if n_threads == 1:
            for s in self.lookup_table:
                all_entries += [ParallelHamiltonianBuilder.do_single_state(s)]

        else:
            with Pool(n_threads) as pool:
                all_entries = pool.map(ParallelHamiltonianBuilder.do_single_state, self.lookup_table)

        # Make a sparse matrix out ot this -although pretty plain, this can handle
        # reasonably sized lists of indices (will do fo now).
        icol, irow, idata = [], [], []
        for line in all_entries:
            if len(line):
                row, col, data = zip(*line)

                # Convert the columns to the proper format.
                for k in range(len(row)):
                    icol.append(col[k])
                    irow.append(self.state_to_index(row[k]))
                    idata.append(data[k])

        if not self.silent:
            self._log("# of nonzero entries: " + str(len(idata)))
        return GaussLatticeHamiltonian(idata, irow, icol, n_fock=self.n_fock)
