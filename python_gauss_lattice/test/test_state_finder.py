""" ----------------------------------------------------------------------------

    test_state_finder.py - LR, May 2020

    Some sanity checks for the state finding algorithm and its utilities.

---------------------------------------------------------------------------- """
from gauss_lattice import GaussLattice
from gauss_lattice.bit_magic import set_bits, sum_occupancies
import numpy as np


def test_winding_links_2D():
    """
    """
    glatt = GaussLattice(L=[2,6])
    e_wx = [0, 4, 8, 12, 16, 20]
    e_wy = [1, 3]
    assert e_wx == sorted(glatt.winding_masks[0])
    assert e_wy == sorted(glatt.winding_masks[1])


def test_winding_links_3D():
    """ Checks whether the correct links are summed for the computation of the
        winding numbers.
    """
    # 2 x 2 x 4
    glatt = GaussLattice(L=[2,2,4])
    e_wx = [0, 6, 12, 18, 24, 30, 36, 42]
    e_wy = [1, 4, 13, 16, 25, 28, 37, 40]
    e_wz = [2, 5, 8, 11]
    assert e_wx == sorted(glatt.winding_masks[0])
    assert e_wy == sorted(glatt.winding_masks[1])
    assert e_wz == sorted(glatt.winding_masks[2])

    # Future.
    # 2 x 4 x 4
    # glatt = GaussLattice(L=[2,4,4])
    # e_wx = [0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84, 90]
    # assert e_wx == sorted(glatt.winding_masks[0])


def check_winding_bins():
    """ Checks whether the correct number of winding bins is assumed.
    """
    glatt = GaussLattice(L=[2,2,4])
    assert glatt.winding_bins.shape == (9, 9, 5)

    glatt = GaussLattice(L=[2,4,2])
    assert glatt.winding_bins.shape == (9, 5, 9)

    glatt = GaussLattice(L=[2,6])
    assert glatt.winding_bins.shape == (13,5)
