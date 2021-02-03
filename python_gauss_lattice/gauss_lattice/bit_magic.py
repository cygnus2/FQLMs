""" ----------------------------------------------------------------------------

    bit_magic.py - LR, May 2020

    Some useful bit routines.

---------------------------------------------------------------------------- """
from numba import jit
import numpy as np

def set_bits(bits):
    if len(bits) != len(set(bits)):
        raise ValueError('Index was provided multiple times.')
    return sum([1<<k for k in bits])

def count_particles(state, nb=64):
    """ Counts particles, i.e., set bits, in an integer. Default: 64 bit, can
        be arbitrary.
    """
    return sum([(state>>k)&1 for k in range(nb)])

# def sum_occupancies(a, b, state):
#     """ Sums all occupancies between index a and b, both inclusive. Order of a
#         and b does not matter, will be from min to max.
#     """
#     return sum([(state>>k)&1 for k in range(min(a,b),max(a,b))])

# @jit(nopython=True)
def sum_occupancies_ordered(a, b, state):
    """ Sums all occupancies between index a and b, both inclusive.
        Assumes a >= b.
    """
    o = 0
    for k in range(b,a+1):
        o += (state>>k)&1
    return o
    # return np.sum([(state>>k)&1 for k in range(b,a)])
