""" ----------------------------------------------------------------------------

    state_finder.py - LR, March 2020

    Finds and stores the Gauss law states for a given lattice configuration.

---------------------------------------------------------------------------- """
from gauss_lattice import GaussLattice
from gauss_lattice.aux import timeit, file_tag, size_tag


def write_winding_sectors(L, wn):
    """ Dirty output function.
    """
    filename='output/winding_sectors_' +size_tag(L) + '.dat'
    if len(L) == 3:
        with open(filename, 'w') as f:
            for i in range(wn.shape[0]):
                for j in range(wn.shape[1]):
                    for k in range(wn.shape[2]):
                        f.write("{:d},{:d},{:d},{:d}\n".format(i,j,k,wn[i,j,k]))


@timeit
def wrap_state_finder(gl, *args, **kwargs):
    """ Wrapper function, purely for timing.
    """
    return gl.find_states(*args, **kwargs)


# Create a GaussLattice with appropriate parameters. If the optional keyword
# `state_file` is provided, *all* states will be stored - either in a hdf5 or
# plain text format, depending of the ending of the filename.
L = [2,2,2]

# glatt = GaussLattice(L=L)
state_file = file_tag(L, filetype='hdf5')
glatt = GaussLattice(L=L, state_file=state_file)

# Constructs the states & times the execution.
wn = wrap_state_finder(glatt)

print(f"Found {wn.sum()} states in total.")
print(f"Found {wn.max()} states in the largest winding sector.")

# Output.
write_winding_sectors(L, wn)
