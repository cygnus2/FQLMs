from .hamiltonian_builder import HamiltonianBuilder
from multiprocessing import Pool
import os, subprocess
import h5py as hdf
import time
import datetime as dt
import numpy as np
from itertools import product
from copy import copy


def _cycle_plaquettes_bare(args):
    # self._log(str(self.level) + " // " + str(state))
    state, plaquettes = args

    states = set()
    for p in plaquettes:
        new_state, _ = apply_u_dagger_bare(state, p, sign=False)
        if not new_state:
            new_state, _ = apply_u_bare(state, p, sign=False)
        if new_state:
            states.add(new_state)
    return states

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


    def level_alert(self, level, number):
        from gauss_lattice.lr_notify import push_message
        push_message(f'Recursion level {level}  with {number} states is completed @{self.host[0]}')

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
        # self._log(str(self.level) + " // " + str(state))
        states = set()
        for p in self.plaquettes:
            new_state, _ = self.apply_u_dagger(state, p, sign=False)
            if not new_state:
                new_state, _ = self.apply_u(state, p, sign=False)
            if new_state:
                states.add(new_state)
        return states


    def _exhaust(self, seed_states, rest, level, pool=None, output_file=None, max_level=10000):
        """ One step in the iteration.
        """
        self.level = level
        self._log(dt.datetime.now().strftime("[%H:%M:%S]")+ "  " + str(level) + " " + str(len(seed_states)))

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
            self.level_alert(level, len(seed_states))

        L = len(seed_states)
        if not L or level>=max_level:
            print(f"Terminated at {level} layers.")
            return rest

        if pool:
            # states = set().union(*pool.map_async(self._cycle_plaquettes, seed_states, chunksize=L//self.n_threads).get())
            # states = set().union(*pool.imap_unordered(self._cycle_plaquettes, seed_states))
            states = set().union(*pool.map(_cycle_plaquettes_bare, product(seed_states, [self.plaquettes])))
        else:
            states = set()
            for s in seed_states:
                states = states | self._cycle_plaquettes(s)

        new_rest = seed_states.union(rest)
        states.difference_update(new_rest)
        return self._exhaust(states, new_rest, level+1, pool, output_file, max_level)


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
        self.n_threads = n_threads
        if n_threads > 1:
            with Pool(n_threads) as p:
                states = self._exhaust(set(seed_states), set(), 0, pool=p, output_file=output_file, max_level=max_level)
        else:
            states = self._exhaust(set(seed_states), set(), 0, pool=None, output_file=output_file, max_level=max_level)

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
