""" ----------------------------------------------------------------------------

    le_hamiltonian_builder.py - LR, October 2020

    Constructs the Hamiltonian only for low-energy gauss-lattice states.

---------------------------------------------------------------------------- """
from .hamiltonian_builder import HamiltonianBuilder
from multiprocessing import Pool
from tqdm import tqdm as tbar
import time

class LowEnergyHamiltonianBuilder(HamiltonianBuilder):
    """ Builds a Hamiltonian with low energy states only.
    """

    def __init__(self, param, logger=None, silent=False):
        """ Here we really do nothing, just set up the object.
        """
        # Call the constructor of the parent but without any states yet - they
        # need to be constructed.
        super().__init__(param, [], logger=logger, silent=True)
        self.silent = silent


    def find_le_states(self, seed_states, iterations):
        """ Takes a number of seed states and applies operations to them to find
            states with single excitations. These will be our states for the
            Hamiltonian.

            Note: very inefficient right now, just a quick test.
        """
        states = self._excite(set(seed_states), iterations)

        self.lookup_table = sorted(list(states))
        self.n_fock = len(self.lookup_table)
        if not self.silent:
            self._log(f"Found {len(self.lookup_table)} states in the low-energy sector.")



    def _excite(self, seed_states, level):
        """ One step in the iteration.
        """
        if not level:
            return seed_states

        states = set()
        for seed_state in seed_states:
            for p in self.plaquettes:

                # First apply the U term.
                new_state, _ = self.apply_u_dagger(seed_state, p)

                # If U term was not successful, try the U^dagger term.
                # (the order could have been switched - there's always only one
                # possibility for overlap to be generated)
                if not new_state:
                    new_state, _ = self.apply_u(seed_state, p)

                if new_state:
                    states.add(new_state)

        return seed_states | self._excite(states, level-1)


    def _cycle_plaquettes(self, state):
        states = set()
        for p in self.plaquettes:
            new_state, _ = self.apply_u_dagger(state, p)
            if not new_state:
                new_state, _ = self.apply_u(state, p)
            if new_state:
                states.add(new_state)
        return states

    def _exhaust(self, seed_states, rest, level, pool):
        """ One step in the iteration.
        """
        if not len(seed_states):
            print(f"Terminated at {level} layers.")
            return rest

        print(level, len(seed_states))
        states = set().union(*pool.map(self._cycle_plaquettes, seed_states))

        new_rest = seed_states.union(rest)
        return self._exhaust(states-new_rest, new_rest, level+1, pool)


    def find_all_states(self, seed_states, n_threads=1):
        """ Finds all states.
        """
        if not self.silent:
            self._log(f"Starting search for states in the low-energy sector on {n_threads} cores.")

        # Execute on multiple processes.
        ts = time.time()
        with Pool(n_threads) as p:
            states = self._exhaust(set(seed_states), set(), 0, p)
        te = time.time()
        if not self.silent:
            self._log('Recursion took  %2.2f ms' % ((te - ts) * 1000))

        self.lookup_table = sorted(list(states))
        self.n_fock = len(self.lookup_table)
        if not self.silent:
            self._log(f"Found {len(self.lookup_table)} states in the low-energy sector.")
        return set(self.lookup_table)
