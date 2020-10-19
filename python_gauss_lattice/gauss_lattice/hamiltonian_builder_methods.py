""" ----------------------------------------------------------------------------

    _hamiltonian_builder_parallel.py - LR, October 2020

    Some routines that need to be adapted for a more efficient parallel
    processing of the loops. Mainly circumvents the copying of large amounts of
    data and avoids pickling the builder object itself (as it is done with
    multiprocessing.Pool).

---------------------------------------------------------------------------- """
from copy import copy


def do_single_state(state, plaquettes):
    """ Constructs matrix-elements for a single Fock state, suitable for parallel
        computing.
    """
    states = []
    for p in plaquettes:
        # First apply the U term.
        new_state, sign = apply_u_dagger(state, p)

        # If U term was not successful, try the U^dagger term.
        # (the order could have been switched - there's always only one
        # possibility for overlap to be generated)
        if not new_state:
            new_state, sign = apply_u(state, p)

        if new_state:
            states.append([state, new_state, sign])

    return states


def apply_u(state, p, sign=True):
    """ Wrapper to apply U
    """
    return _apply_plaquette_operator(state, p, mask=[False, False, True, True], sign=sign)

def apply_u_dagger(state, p, sign=True):
    """ Wrapper to apply U+
    """
    return _apply_plaquette_operator(state, p, mask=[True, True, False, False], sign=sign)

def _apply_plaquette_operator(state, p, mask, sign):
    """ Applies the U operator

            c1+ c2+ c3 c4

        or the U+ term

            c1 c2 c3+ c4+

        to a given plaquette in a given state (link configuration starting
        with x-mu link and going counter-clockwise).
    """
    new_state = copy(state)
    for k in range(4):
        m = 1 << p[k]
        if bool(new_state & m) == mask[k]:
            new_state = copy(new_state^m)
        else:
            return 0, 0
    return new_state, 1
