""" ----------------------------------------------------------------------------

    test_defs.py - LR, May 2020

    Some useful definitions for the unittests.

---------------------------------------------------------------------------- """


def _set_bits(bits):
    if len(bits) != len(set(bits)):
        raise ValueError('Index was provided multiple times.')
    return sum([1<<k for k in bits])
