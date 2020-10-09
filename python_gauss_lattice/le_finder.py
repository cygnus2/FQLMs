""" ----------------------------------------------------------------------------

    diagonalizer.py - LR, May 2020

    Run script for diagonalizing the Hamiltonian for the gauss lattice.

---------------------------------------------------------------------------- """
from gauss_lattice import LowEnergyHamiltonianBuilder
from gauss_lattice.aux import load_config, timeit, read_winding_sector, size_tag
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

# These are lattices that are maximally flippable for each lattice size. It does not
# matter where we start, as long as we're in the correct winding sector (and we
# have flippable plaquettes at all).
base_lattices = {
    (2,2,2) : [3816540, 3872106, 5421780, 5678001, 7542990, 7743645,
                9033570, 9234225, 11099214, 11355435, 12905109, 12960675],
    (2,2,4) : [64030919769180, 64963162609002, 90962379586260, 95261054903217, 126550380058830, 129916812535965,
                151558164174690, 154924596651825, 186213921807438, 190512597124395, 216511814101653, 217444056941475],
    (2,2,6) : [1074260571646206819420, 1089901011134353970538,  1526095490192680073940, 1598215294499136381873,
                2123163061129091160270, 2179642425947400317085, 2542724056922244896610, 2599203421740554053425,
                3124151188370508831822, 3196270992676965139755, 3632465471735291243157, 3648105911223438394275],
}

@timeit(logger=None)
def le_finding(builder, *args, **kwargs):
    """ This method just exists to be able to time the routine efficiently.
    """
    return builder.find_all_states(*args, **kwargs)

# Recursively find all states for the
states = le_finding(
    builder,
    base_lattices[tuple(param['L'])],
    n_threads=param["n_threads"],
    output_file=param['working_directory'] + "/le_states_" + size_tag(param['L']) + '.hdf5',
    max_level=param['maximum_excitation_level']
)
if param["L"] == [2,2,2]:
    all_states, _ = read_winding_sector(param['L'], (0,0,0), basedir=param['working_directory'])
    print("States that are not found:")
    print(set(all_states) - states)
