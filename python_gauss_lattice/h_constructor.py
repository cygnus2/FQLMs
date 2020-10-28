""" ----------------------------------------------------------------------------

    h_constructor.py - LR, October 2020

---------------------------------------------------------------------------- """
from gauss_lattice import LowEnergyHamiltonianBuilder, ParallelHamiltonianBuilder, HamiltonianBuilder
from gauss_lattice.aux import size_tag, timeit, read_winding_sector, read_all_states, load_config
import numpy as np
import argparse, logging
import h5py as hdf


@timeit(logger=None)
def hamiltonian_construction(builder, *args, **kwargs):
    """ This method just exists to be able to time the routine efficiently.
    """
    # return builder.construct()
    return builder.construct(*args, **kwargs)

# ------------------------------------------------------------------------------
# Input handling.
parser = argparse.ArgumentParser(description="Python gauss lattice diagonalizer (single lambda).")
parser.add_argument('-i', metavar='', type=str, default=None, help='YAML style input file.')
args = parser.parse_args()

# Get parameters.
param = load_config(args.i)

# Set up a logger with a handler for the terminal output.
logger = logging.getLogger('diagonalization logger')
logger.addHandler(logging.StreamHandler())
logger.addHandler(logging.FileHandler(param['working_directory'] + "/" + param['logfile'], mode='w'))
logger.setLevel(logging.DEBUG)

# ------------------------------------------------------------------------------
# Hamiltonian setup.
if param.get('winding_sector'):
    states, ws = read_winding_sector(param['L'], param['winding_sector'], basedir=param['working_directory'])
else:
    states = read_all_states(param['L'], basedir=param['working_directory'])
    ws = None

builder = ParallelHamiltonianBuilder(param, states, logger=logger)
# builder.read_le_states(states)
ham = hamiltonian_construction(builder, n_threads=param['n_threads'])

# hamiltonian_file = param['working_directory'] + '/SEQUENTIAL_hamiltonian_' + size_tag(param['L']) + '.hdf5'
# If specified, store the Hamiltonian for later use.
# if param['store_hamiltonian']:
#     with hdf.File(hamiltonian_file, 'a' if i else 'w') as f:
#         ds = f.create_dataset(ws, data= np.array([ham.col, ham.row, ham.data]))
#         ds.attrs['n_fock'] = len(states)
