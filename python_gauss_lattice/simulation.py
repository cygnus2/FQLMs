""" ----------------------------------------------------------------------------

    simulation.py - LR, October 2020

    Provides a simulation object for convenience when writing similar yet
    slightly different run scripts.

---------------------------------------------------------------------------- """
import argparse, logging, os
import h5py as hdf
import numpy as np
from gauss_lattice.aux import timestamp, load_config, timeit, file_tag, size_tag, read_winding_sector, read_all_states
from gauss_lattice import GaussLatticeHamiltonian, HamiltonianBuilder


@timeit(logger=None)
def hamiltonian_construction(builder, *args, **kwargs):
    """ This method just exists to be able to time the routine efficiently.
    """
    return builder.construct(*args, **kwargs)

@timeit(logger=None)
def hamiltonian_diagonalization(ham, **kwargs):
    """ This method just exists to be able to time the routine efficiently.
    """
    return ham.diagonalize(**kwargs)


class GLSimulation(object):
    """ Convenience container for repetitive tasks.
    """
    def __init__(self):
        """ Set up the most common structure, gets information directly from the
            parsing of the command line arguments.
        """
        # Translate the arguments.
        parser = argparse.ArgumentParser(description="ED Simulation for gauss lattice.")
        parser.add_argument('-i', metavar='', type=str, default=None, help='YAML style input file.')
        parser.add_argument('-notify', action='store_true', help='Toggles whether a notification should be sent when done. Should only be used by Lukas (sorry).')
        args = parser.parse_args()

        # This sets the parameters fof the calculation (everything else is fixed).
        self.param = load_config(args.i)
        self.notify = args.notify

        # Make output destination.
        self.working_directory = self.param.get('working_directory', 'output/')
        os.makedirs(self.working_directory, exist_ok=True)

        # Set up a logger with a handler for the terminal output.
        self.logger = logging.getLogger('parameter run logger')
        self.logger.addHandler(logging.StreamHandler())
        self.logger.addHandler(logging.FileHandler(self.working_directory + '/' + self.param['logfile'], mode='w'))
        self.logger.setLevel(logging.DEBUG)

        self.log(f'Working in directory {self.working_directory}')

        # Set some further flags.
        self.compute_eigenstates = self.param.get('compute_eigenstates', False)


    def log(self, msg):
        self.logger.info(timestamp() + ' ' + msg)


    def _get_state_file(self, default=None):
        """ Getter for the file where states are stored/loaded.
        """
        state_file = self.param.get("state_file", default)
        if not state_file:
            state_file = self.working_directory+file_tag(self.param['L'], filetype='hdf5').replace("winding_", "le_")
        return state_file


    def _get_hamiltonian_file(self, default=None):
        """ Getter for the file where the hamiltonian is stored/loaded.
        """
        hamiltonian_file = self.param.get('hamiltonian_file', default)
        if not hamiltonian_file:
            hamiltonian_file = self.working_directory + '/hamiltonian_' + size_tag(self.param['L']) + '.hdf5'
        return hamiltonian_file


    def _get_result_file(self, prefactor='results', default=None):
        """ File where a spectrum is stored.
        """
        result_file = self.param.get('result_file', default)
        if not result_file:
            result_file = (
                self.working_directory +
                '/'+prefactor+'_'+
                self.param['gauge_particles'] + '_' +
                ('' if self.ws is None else self.ws+'_') +
                size_tag(self.param['L']) +
                '.hdf5'
            )
        return result_file


    def read_states(self, file=None):
        """ Read the GLS from file.
        """
        state_file = self._get_state_file(default=file)
        try:
            if self.param.get('winding_sector'):
                states, self.ws = read_winding_sector(
                    self.param['L'],
                    self.param['winding_sector'],
                    filename=state_File
                )
            else:
                states = read_all_states(
                    self.param['L'],
                    filename=state_file
                )
                self.ws = 'all-ws'
            self.log(f'Read Fock states from {state_file}')
            return states

        except (OSError):
            self.log(f'Could not find file {state_file}')
            return None


    def read_hamiltonian(self, ham_name, file=None):
        """ Read the Hamiltonian from file.
        """
        ham_file = self._get_hamiltonian_file(default=file)
        try:
            with hdf.File(ham_file, 'r') as f:
                mat = f[ham_name][...]
                ham = GaussLatticeHamiltonian(mat[2,:], mat[1,:], mat[0,:], f[ham_name].attrs['n_fock'])
            self.log(f'Read Hamiltonian from {ham_file}')
            return ham

        except (OSError):
            self.log(f'Could not find file {ham_file}')
            return None

        except (KeyError, ValueError):
            self.log(f'Could not find object {ham_name} in file {ham_file}')
            return None


    def construct_hamiltonian(self, states, ham_type=HamiltonianBuilder, **kwargs):
        """ Takes a list of states and constructs the Hamiltonian.
        """
        builder = ham_type(
            self.param,
            states=states,
            logger=self.logger,
            **kwargs
        )
        self.store_ham = self.param.get('store_hamiltonian', False)
        return hamiltonian_construction(builder, self.param.get('n_threads', 1))


    def diagonalize_hamiltonian(self, ham, lam=None):
        """ Diagonalizes the Hamiltonian with given parameters.
        """
        return hamiltonian_diagonalization(
            ham,
            full_diag = self.param.get('full_diag'),
            J = self.param['J'],
            lam = lam if lam else self.param['lambda'],
            gauge_particles = self.param['gauge_particles'],
            n_eigenvalues = max(1, self.param['n_eigenvalues']),
            which = self.param.get('ev_type', 'SA'),
            compute_eigenstates = self.compute_eigenstates
        )


    def store_hamiltonian(self, ham, label, file=None, grp_name=None, attrs={}):
        """ Stores the Hamiltonian.
        """
        ham_file = file if file else self._get_hamiltonian_file()
        self.log(f'Storing Hamiltonian in {ham_file}')

        self.store_data(
            np.array([ham.col, ham.row, ham.data]),
            label,
            grp_name=grp_name,
            attrs=attrs,
            file=ham_file
        )


    def store_data(self, data, ds_name, grp_name=None, attrs={}, file=None):
        """ Stores the results in standardized fashion.
        """
        filename = file if file else self._get_result_file(default=file)
        with hdf.File(filename, 'a') as f:
            ds = None

            # Store in group or directly in dataset.
            if grp_name:
                grp = f[grp_name] if grp_name in f else f.create_group(grp_name)
                if ds_name in grp:
                    del grp[ds_name]
                ds = grp.create_dataset(ds_name, data=data)
            else:
                if ds_name in f:
                    del f[ds_name]
                ds = f.create_dataset(ds_name, data=data)

            # Finally, store attributes.
            for k, v in attrs.items():
                ds.attrs[k] = v


    def push_notification(self, msg):
        """ If desired, a notification will be sent to Lukas' phone.
        """
        if self.notify:
            from gauss_lattice.lr_notify import push_message
            import os, subprocess
            host = subprocess.check_output(['hostname']).strip().decode('UTF-8'),
            push_message(f'[{host[0]}] ' + msg)
