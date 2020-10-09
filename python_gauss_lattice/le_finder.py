""" ----------------------------------------------------------------------------

    diagonalizer.py - LR, May 2020

    Run script for diagonalizing the Hamiltonian for the gauss lattice.

---------------------------------------------------------------------------- """
from gauss_lattice import LowEnergyHamiltonianBuilder
from gauss_lattice.aux import load_config, timeit, read_winding_sector
import argparse, logging

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
builder = LowEnergyHamiltonianBuilder(param, logger=logger)
base_lattices = {
    (2,2,2) : [3816540],
    #
    # (2,2,2) : [3816540, 3872106, 5421780, 5678001, 7542990, 7743645,
    #             9033570, 9234225, 11099214, 11355435, 12905109, 12960675],
    (2,2,4) : [64030919769180, 64963162609002, 90962379586260, 95261054903217,
                126550380058830, 129916812535965, 151558164174690, 154924596651825,
                186213921807438, 190512597124395, 216511814101653, 217444056941475]
 }



@timeit(logger=None)
def le_finding(builder, *args, **kwargs):
    """ This method just exists to be able to time the routine efficiently.
    """
    return builder.find_all_states(*args, **kwargs)

# Recursively find all states for the
states = le_finding(builder, base_lattices[tuple(param['L'])], n_threads=param["n_threads"])
if param["L"] == [2,2,2]:
    all_states, _ = read_winding_sector(param['L'], (0,0,0), basedir=param['working_directory'])
    print("States that are not found:")
    print(set(all_states) - states)
