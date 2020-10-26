""" ----------------------------------------------------------------------------

    _hamiltonian_builder_parallel.py - LR, October 2020

    Some routines that need to be adapted for a more efficient parallel
    processing of the loops. Mainly circumvents the copying of large amounts of
    data and avoids pickling the builder object itself (as it is done with
    multiprocessing.Pool).

---------------------------------------------------------------------------- """
from copy import copy
from .bit_magic import sum_occupancies_ordered
import numpy as np


def do_single_state(args, sign=True, set_collection=False):
    """ Constructs matrix-elements for a single Fock state, suitable for parallel
        computing.
    """
    state, plaquettes = args

    states = set() if set_collection else []
    for p in plaquettes:
        # First apply the U term.
        new_state, s = apply_u_dagger(state, p, sign=sign)

        # If U term was not successful, try the U^dagger term.
        # (the order could have been switched - there's always only one
        # possibility for overlap to be generated)
        if not new_state:
            new_state, s = apply_u(state, p, sign=sign)

        if new_state:
            if set_collection:
                states.add(new_state)
            else:
                states.append([state, new_state, s])

    return states



def cycle_plaquettes(args):
    return do_single_state(args, sign=False, set_collection=True)

def apply_u(state, p, sign=True):
    """ Wrapper to apply U
    """
    return _apply_plaquette_operator(state, p, [False, False, True, True], sign)

def apply_u_dagger(state, p, sign=True):
    """ Wrapper to apply U+
    """
    return _apply_plaquette_operator(state, p, [True, True, False, False], sign)


def _apply_plaquette_operator(state, p, mask, sign):
    """ Applies the U operator

            c1+ c2+ c3 c4

        or the U+ term

            c1 c2 c3+ c4+

        to a given plaquette in a given state (link configuration starting
        with x-mu link and going counter-clockwise).
    """
    n = 0
    if sign:
        a = max(p[:-1])
    new_state = copy(state)
    for k in range(4):
        m = 1 << p[k]
        if bool(new_state & m) == mask[k]:
            if sign:
                n += sum_occupancies_ordered(a, p[k], new_state)
            new_state = copy(new_state^m)
        else:
            return 0, 0

    return new_state, (-1)**n
