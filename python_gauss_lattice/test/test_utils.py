""" ----------------------------------------------------------------------------

    test_utils.py - LR, May 2020

    Ensure the utilities work well.

---------------------------------------------------------------------------- """
import numpy as np
import pytest
from gauss_lattice.bit_magic import set_bits


def test_state_gen():
    """ Quick test fo state generator.
    """
    bits = set(np.random.randint(32, size=np.random.randint(5, 10)))
    assert set_bits(bits) == sum([2**k for k in bits])

    with pytest.raises(ValueError):
        state = set_bits([1,1,2,3])
