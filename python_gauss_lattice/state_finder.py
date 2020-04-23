from gl_aux import timeit
from gauss_lattice import GaussLattice

@timeit
def wrap_state_finder(gl, *args, **kwargs):
    return gl.find_states(*args, **kwargs)


# Create a GaussLattice with appropriate parameters.
L = [2,8]
# glatt = GaussLattice(L=L)
glatt = GaussLattice(L=L, state_file='test.dat', buffer_length=1000)

# Constructs the states & times the execution.
wn = wrap_state_finder(glatt)
print(wn)
print(f"Found {wn.sum()} states in total.")
print(f"Found {wn.max()} states in the larges winding sector.")


if len(L) == 3:
    with open('output/states_full.txt', 'w') as f:
        for i in range(wn.shape[0]):
            for j in range(wn.shape[1]):
                for k in range(wn.shape[2]):
                    f.write("{:d},{:d},{:d},{:d}\n".format(i,j,k,wn[i,j,k]))
