""" ----------------------------------------------------------------------------

    bit_magic.py - LR, May 2020

    Some useful bit routines.

---------------------------------------------------------------------------- """
def set_bits(bits):
    if len(bits) != len(set(bits)):
        raise ValueError('Index was provided multiple times.')
    return sum([1<<k for k in bits])

def count_particles(state, nb=64):
    """ Counts particles, i.e., set bits, in an integer. Default: 64 bit, can
        be arbitrary.
    """
    return sum([(state>>k)&1 for k in range(nb)])

def sum_occupancies(a, b, state):
    """ Sums all occupancies between index a and b, both inclusive. Order of a
        and be does not matter, will be from min to max.
    """
    return sum([(state>>k)&1 for k in range(min(a,b),max(a,b))])
