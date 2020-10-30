""" ----------------------------------------------------------------------------

    gl_simulation.py - LR, October 2020

    Provides a simulation object for convenience when writing similar yet
    slightly different run scripts.

---------------------------------------------------------------------------- """
import argparse, logging, os, yaml, sys
import subprocess
import h5py as hdf
import numpy as np
from copy import copy

sys.path.append('../')
from gauss_lattice.aux import timeit, timestamp, full_timestamp
from gauss_lattice import GaussLatticeHamiltonian, HamiltonianBuilder, GaussLattice


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
        self.param = self._load_config(args.i)
        self.notify = args.notify

        # Make output destination.
        self.working_directory = self.param.get('working_directory', 'output/')
        os.makedirs(self.working_directory, exist_ok=True)

        # Set up a logger with a handler for the terminal output.
        self.logger = logging.getLogger('parameter run logger')
        self.logger.addHandler(logging.StreamHandler())
        self.logger.addHandler(logging.FileHandler(self.working_directory + '/' + self.param['logfile'], mode='w'))
        self.logger.setLevel(logging.DEBUG)

        self.host = subprocess.check_output(['hostname']).strip().decode('UTF-8')
        self.log(f'Working in directory {self.working_directory} at host {self.host}')


        # Set the winding sector (mainly for output and correct state generation).
        if 'winding_sector' in self.param:
            # The winding sectors are labelled differently in the HDF5 file (for
            # convenience reasons) so we have to relabel them here.
            ws_shifted = GLSimulation._winding_shift(self.param['L'], self.param['winding_sector'])
            self.ws = GLSimulation._winding_tag(ws_shifted)
        else:
            self.ws = 'all-ws'

        # Set some further flags.
        self.compute_eigenstates = self.param.get('compute_eigenstates', False)


    # --------------------------------------------------------------------------
    # Actual calculation routines.

    def find_states(self, *args, file=None, **kwargs):
        """ Find the GL states with recursive algorithm (slower than with LE),
            only recommended for smaller runs (maximally 2x2x4).

            The other use-case is to provide a file with the states.
        """
        self.log('Could not find stored states, constructing Fock state list from scratch.')
        state_file = self._get_state_file(default=file)
        glatt = GaussLattice(
            self.param['L'],
            state_file=state_file,
            filetype='hdf5',
            basedir='/' if self.working_directory[0]=='/' else './'
        )
        glatt.find_states()

        # Read from file again, cumbersome but this is the legacy structure.
        return self.read_states(*args, file=file, **kwargs)

    def construct_hamiltonian(self, states, builder_type=HamiltonianBuilder, **kwargs):
        """ Takes a list of states and constructs the Hamiltonian.
        """
        builder = builder_type(
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


    def run_lambda_loop(self, lambdas, ham):
        """ Wrapper to run multiple values of lambdas back-to-back.
        """
        for i, l in enumerate(lambdas):
            self.log('[{:d} / {:d}] diagonalizing Hamiltonian for lambda={:.4f}'.format(i+1, len(lambdas), l))

            # Diagonalization.
            results = self.diagonalize_hamiltonian(ham, lam=l)
            attrs = {
                'lambda' : l
            }
            if not self.compute_eigenstates:
                self.store_data(results, ds_name='spectrum_lam_{:6f}'.format(l), attrs=attrs)
            else:
                self.store_data(results[0], ds_name='spectrum_lam_{:6f}'.format(l), attrs=attrs)
                self.store_data(results[1], ds_name='eigenstates_lam_{:6f}'.format(l), attrs=attrs)


    # --------------------------------------------------------------------------
    # I/O stuff.

    def _load_config(self, filename):
        """ Reads a config from a yaml file, with appropriate chekcs.
        """
        if filename is not None:
            # Loads parameters from file.
            with open(filename, 'r') as f:
                try:
                    return yaml.safe_load(f)
                except yaml.YAMLError as exc:
                    print(exc)
                    raise yaml.YAMLError()
        else:
            sys.exit('fatal: no input file specified')


    def log(self, msg):
        self.logger.info(timestamp() + ' ' + msg)


    def _get_state_file(self, prefactor='winding_states', default=None):
        """ Getter for the file where states are stored/loaded.
        """
        state_file = self.param.get("state_file", default)
        if not state_file:
            state_file = (
                self.working_directory +
                '/' + prefactor + '_' +
                GLSimulation._size_tag(self.param['L']) +
                '.hdf5'
            )
        return state_file

    def _get_hamiltonian_file(self, prefactor='hamiltonian', default=None):
        """ Getter for the file where the hamiltonian is stored/loaded.
        """
        hamiltonian_file = self.param.get('hamiltonian_file', default)
        if not hamiltonian_file:
            hamiltonian_file = (
            self.working_directory +
            '/' + prefactor + '_' +
            GLSimulation._size_tag(self.param['L']) +
            '.hdf5'
        )
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
                GLSimulation._size_tag(self.param['L']) +
                '.hdf5'
            )
        return result_file


    def read_states(self, file=None, merged=True):
        """ Read the GLS from file.
        """
        state_file = self._get_state_file(default=file)
        try:
            states = []
            if self.param.get('winding_sector'):
                with hdf.File(state_file, 'r') as f:
                    states = f[self.ws][...]

            else:
                with hdf.File(state_file, 'r') as f:
                    for ws in f:
                        if merged:
                            states += list(f[self.ws][...])
                        else:
                            states.append([ws, list(f[self.ws][...])])
            self.log(f'Read Fock states from {state_file}')
            return states

        except (OSError):
            self.log(f'Could not find file {state_file}')
            return []


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

        all_attrs = copy(attrs)
        all_attrs['time'] = full_timestamp()
        all_attrs['host'] = self.host

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
            for k, v in all_attrs.items():
                ds.attrs[k] = v

    # --------------------------------------------------------------------------
    # Auxiliary stuff.

    def push_notification(self, msg):
        """ If desired, a notification will be sent to Lukas' phone.
        """
        if self.notify:
            from gauss_lattice.lr_notify import push_message
            push_message(f'[{self.host}] ' + msg)


    @staticmethod
    def _winding_tag(ws, labels=['x', 'y', 'z'], shift=None):
        """ Returns the naming convention of the winding datasets.

            The shift is a lattice configuration L = [Lx,Ly,Lz], such
            that the HDF5 datasets may be resolved.
        """
        if shift is not None:
            ws = _winding_shift(shift, ws)
        wtag = ''
        for k in range(len(ws)):
            wtag += 'w{:s}_{:d}-'.format(labels[k], ws[k])
        return wtag[:-1]


    @staticmethod
    def _winding_shift(L, ws):
        """ Maps between the representation of winding numbers.
        """
        if len(L) == 2:
            shift = np.array(L[::-1]) // 2
        elif len(L) == 3:
            shift = np.array([
                L[1]*L[2] // 2,
                L[0]*L[2] // 2,
                L[0]*L[1] // 2
            ])
        else:
            raise NotImplementedError('Dimension not implemented!')
        return ws + shift

    @staticmethod
    def _size_tag(L):
        stag = ''
        for k in range(len(L)):
            stag += '{:d}x'.format(L[k])
        return stag[:-1]
