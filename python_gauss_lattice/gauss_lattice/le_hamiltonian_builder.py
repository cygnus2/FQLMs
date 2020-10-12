""" ----------------------------------------------------------------------------

    le_hamiltonian_builder.py - LR, October 2020

    Constructs the Hamiltonian only for low-energy gauss-lattice states.

---------------------------------------------------------------------------- """
from .hamiltonian_builder import HamiltonianBuilder
from multiprocessing import Pool
import os, subprocess
import h5py as hdf
import time
import numpy as np

class LowEnergyHamiltonianBuilder(HamiltonianBuilder):
    """ Builds a Hamiltonian with low energy states only.
    """

    def __init__(self, param, logger=None, silent=False, notify_level=0):
        """ Here we really do nothing, just set up the object.
        """
        # Call the constructor of the parent but without any states yet - they
        # need to be constructed.
        super().__init__(param, [], logger=logger, silent=True)
        self.silent = silent

        self.notify_level = notify_level
        if notify_level:
            self.host = subprocess.check_output(['hostname']).strip().decode('UTF-8'),


    def level_alert(self, level):
        from gauss_lattice.lr_notify import push_message
        push_message(f'Recursion level {level}  @{self.host[0]} is exported!')

    def _split(self, data, bit_shift=63):
        """ Converts to two integer lists to be able to store it in HDF5 (which
            can maximally do 64 bit numbers).
        """
        div = 2**bit_shift
        split_data = np.zeros((len(data),2), dtype=np.int64)
        for i, x in enumerate(data):
            split_data[i,:] = x >> bit_shift, x%div
        return split_data


    def _cycle_plaquettes(self, state):
        states = set()
        for p in self.plaquettes:
            new_state, _ = self.apply_u_dagger(state, p)
            if not new_state:
                new_state, _ = self.apply_u(state, p)
            if new_state:
                states.add(new_state)
        return states


    def _exhaust(self, seed_states, rest, level, pool, output_file=None,max_level=10000):
        """ One step in the iteration.
        """
        self._log(str(level) + " " + str(len(seed_states)))

        # Write states.
        if output_file:
            with hdf.File(output_file, 'r+') as f:
                if self.big_int:
                    data = self._split(list(seed_states))
                else:
                    data = list(seed_states)
                f.create_dataset('states_lv_{:d}'.format(level), data=data)

        # Notify for testing purposes.
        if self.notify_level and (level>=self.notify_level):
            self.level_alert(level)

        if not len(seed_states) or level>=max_level:
            print(f"Terminated at {level} layers.")
            return rest

        states = set().union(*pool.map(self._cycle_plaquettes, seed_states))

        new_rest = seed_states.union(rest)
        return self._exhaust(states-new_rest, new_rest, level+1, pool, output_file,max_level)


    def find_all_states(self, seed_states, n_threads=1, output_file=None, max_level=10000):
        """ Finds all states.
        """
        if not self.silent:
            self._log(f"Starting search for states in the low-energy sector on {n_threads} cores.")

        # Prepare
        if output_file:
            with hdf.File(output_file, 'w') as f:
                f.attrs['n_threads'] = n_threads
            self._log(f"Writing states to {output_file}")

        # Execute on multiple processes.
        ts = time.time()
        with Pool(n_threads) as p:
            states = self._exhaust(set(seed_states), set(), 0, p, output_file,max_level=max_level)
        te = time.time()
        if not self.silent:
            self._log('Recursion took  %2.2f ms' % ((te - ts) * 1000))

        self.lookup_table = sorted(list(states))
        self.n_fock = len(self.lookup_table)
        if not self.silent:
            self._log(f"Found {len(self.lookup_table)} states in the low-energy sector.")
        return set(self.lookup_table)


    def _combine(self, data, bit_shift=63):
        print(len(data))
        combined_data = [2**70]*len(data)
        for i, (x, y) in enumerate(data):
            combined_data[i] = (int(x)<<bit_shift) + int(y)
        return sorted(combined_data)

    def read_le_states(self, states):
        """ Takes a list of states and translates them. If the integer size is
            exceeded, the states will be reconstructed from two integers.
        """
        if self.big_int:
            self.lookup_table = self._combine(states)
        else:
            self.lookup_table = sorted(states)

        self.n_fock = len(self.lookup_table)
        if not self.silent:
            self._log(f"Read {self.n_fock} states.")
