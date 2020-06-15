""" ----------------------------------------------------------------------------

    state_finder.py - LR, March 2020

    Finds and stores the Gauss law states for a given lattice configuration.

---------------------------------------------------------------------------- """
import argparse
from gauss_lattice import GaussLattice
from gauss_lattice.aux import timeit, file_tag, size_tag, winding_sectors, load_config


def write_winding_sectors(L, wn, basedir):
    """ Dirty output function.
    """
    filename=basedir+'winding_sectors_' +size_tag(L) + '.dat'
    with open(filename, 'w') as f:
        for ws in winding_sectors(L):
            f.write(
                ("{:d},"*len(L)).format(*ws) +
                "{:d}\n".format(wn[ws])
            )

@timeit(logger=None)
def wrap_state_finder(gl, *args, **kwargs):
    """ Wrapper function, purely for timing.
    """
    return gl.find_states(*args, **kwargs)


parser = argparse.ArgumentParser(description="Python gauss lattice state finder.")
parser.add_argument('-i', metavar='', type=str, default=None, help='YAML style input file (must hold \'L\' and \'working_directory\').')
parser.add_argument('-notify', action='store_true', help='Toggles whether a notification should be sent when done. Should only be used by Lukas (sorry).')
args = parser.parse_args()

# This sets the parameters fof the calculation (everything else is fixed).
param = load_config(args.i)

# Create a GaussLattice with appropriate parameters & stores the states..
state_file = file_tag(param['L'], filetype='hdf5')
glatt = GaussLattice(L=param['L'], state_file=state_file, basedir=param['working_directory'])

# Constructs the states & times the execution.
wn = wrap_state_finder(glatt)

print(f"Found {wn.sum()} states in total.")
print(f"Found {wn.max()} states in the largest winding sector.")

# Output.
write_winding_sectors(param['L'], wn, param['working_directory'])
