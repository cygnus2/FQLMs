from .hamiltonian_builder import HamiltonianBuilder
from .hamiltonian_builder_methods import cycle_plaquettes
from .aux_stuff import timestamp
from multiprocessing import Pool
import os, subprocess
import h5py as hdf
import time
import datetime as dt
import numpy as np
from itertools import product
from copy import copy


class LowEnergyStateFinder(HamiltonianBuilder):
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


    def _exhaust(self, seed_states, rest, level, pool=None, output_file=None, max_level=10000):
        """ One step in the iteration.
        """
        self.level = level
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
            self.level_alert(level, len(seed_states))

        L = len(seed_states)
        if not L or level>=max_level:
            self._log(f"Terminated at {level} layers.")
            return seed_states.union(rest)

        if pool:
            states = set().union(*pool.map(cycle_plaquettes, product(seed_states, [self.plaquettes])))
        else:
            # states = set()
            # for s in seed_states:
            #     states = states | cycle_plaquettes((s,self.plaquettes))
            setlist = [cycle_plaquettes((s,self.plaquettes)) for s in seed_states]
            states = set().union(*setlist)

        new_rest = seed_states.union(rest)
        states.difference_update(new_rest)
        # states = states.difference(new_rest)
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
