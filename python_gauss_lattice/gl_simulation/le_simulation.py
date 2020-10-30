""" ----------------------------------------------------------------------------

    le_simulation.py - LR, October 2020

    Low-energy version of the simulation object.

---------------------------------------------------------------------------- """
import sys
import numpy as np
import h5py as hdf
from .gl_simulation import GLSimulation

sys.path.append('../')
from gauss_lattice import LowEnergyStateFinder


class LowEnergyGLSimulation(GLSimulation):
    """ Extenstion of the Simulation class for the specific purpose of low
        energy runs.
    """
    def __init__(self, *args, **kwargs):
        GLSimulation.__init__(self, *args, **kwargs)
        self.ws = 'all-ws'

    @staticmethod
    def _combine(data, bit_shift=63):
        """ Converter method required to map back two integers with maximum
            64 bit (like HDF5 has it) to an unlimited integer of Python. This
            is needed for lattices with more than 63 spins.
        """
        combined_data = [2**70]*len(data)
        for i, (x, y) in enumerate(data):
            combined_data[i] = (int(x)<<bit_shift) + int(y)
        return sorted(combined_data)


    def _get_hamiltonian_file(self, prefactor='LE_hamiltonian', default=None):
        return GLSimulation._get_hamiltonian_file(self, prefactor=prefactor, default=default)

    def _get_result_file(self, prefactor='LE_results', default=None):
        return GLSimulation._get_result_file(self, prefactor=prefactor, default=default)

    def _get_state_file(self, prefactor='le_states', default=None):
        return GLSimulation._get_state_file(self, prefactor=prefactor, default=default)


    def read_le_states(self, max_level, combine=True, file=None):
        """ Reads all states up to a certain level.
        """
        state_file = self._get_state_file(default=file)
        try :
            states = {}
            with hdf.File(state_file, 'r') as f:
                levels = []
                for g in f:
                    l = int(g.split("_")[-1])
                    levels.append(l)
                    if l <= max_level:
                        s = f[g][...]
                        if len(s.shape) == 2:
                            states[l] = LowEnergyGLSimulation._combine(s)
                        else:
                            states[l] = sorted(map(int, states))

                if sorted(levels)[-1] < max_level:
                    self.log(f"Warning: maximal level not reachable from list! Taking maximum of {l}.")
                    new_max = l
                else:
                    new_max = max_level


            self.log(f'Read Fock states from {state_file}')
            if combine:
                return np.concatenate([states[k] for k in sorted(states.keys()) if len(states[k])]), new_max
            return states, new_max

        except OSError:
            self.log(f'Could not find file {state_file}')
            return []


    def find_le_states(self, base_lattices):
        """ Finds the low energy states starting from the base states provided.
        """
        self.log('Could not find stored states, constructing low-energy Fock state list from scratch.')

        lesf = LowEnergyStateFinder(self.param, self.logger)
        states = lesf.find_all_states(
            base_lattices,
            n_threads=self.param.get('n_threads', 1),
            max_level=self.param.get('maximum_excitation_level', 10000)
        )
        return states


    def store_le_states(self):
        raise NotImplementedError('On-the-fly state storage is not implemented yet!')
