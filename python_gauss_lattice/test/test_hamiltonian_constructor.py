""" ----------------------------------------------------------------------------

    test_hamiltonian_constructor.py - LR, May 2020

    Some sanity checks for the Hamiltonian constructor.

---------------------------------------------------------------------------- """
from gauss_lattice import HamiltonianBuilder
from test_defs import _set_bits
import numpy as np


def test_hamiltonian_builder_index_shift_2D():
    """ Check if PBC are enforced correctly in the Hamiltonian construction (in
        fact, in the construction of the plaquette list).
    """
    # 2D check.
    param = {
        'L' : [2,4]
    }
    builder = HamiltonianBuilder(param, states=[])

    # Set up test case.
    expected = np.array([
        [1, 2],
        [0, 3],
        [3, 4],
        [2, 5],

        [5, 6],
        [4, 7],
        [7, 0],
        [6, 1],
    ])
    for n in range(builder.S[-1]):
        for d in range(len(param['L'])):
            assert expected[n,d] == builder.shift_index(n, d)



def test_hamiltonian_builder_index_shift_3D():
    """ Check if PBC are enforced correctly in the Hamiltonian construction (in
        fact, in the construction of the plaquette list).
    """
    param = {
        'L' : [2,2,4]
    }
    builder = HamiltonianBuilder(param, states=[])

    # Set up test case.
    expected = np.array([
        [1, 2, 4],
        [0, 3, 5],
        [3, 0, 6],
        [2, 1, 7],

        [5, 6, 8],
        [4, 7, 9],
        [7, 4, 10],
        [6, 5, 11],

        [9, 10, 12],
        [8, 11, 13],
        [11, 8, 14],
        [10, 9, 15],

        [13, 14, 0],
        [12, 15, 1],
        [15, 12, 2],
        [14, 13, 3],
    ])
    for n in range(builder.S[-1]):
        print(n)
        for d in range(len(param['L'])):
            assert expected[n,d] == builder.shift_index(n, d)



def test_plaquette_list_2D():
    """ Check if the correct list of plaquettes is found in 2D.
    """
    param = {
        'L' : [2,4]
    }
    builder = HamiltonianBuilder(param, states=[])

    # Set up test case.
    expected = [
        [0, 3, 4, 1],
        [2, 1, 6, 3],
        [4, 7, 8, 5],
        [6, 5, 10, 7],

        [8, 11, 12, 9],
        [10, 9, 14, 11],
        [12, 15, 0, 13],
        [14, 13, 2, 15],
    ]

    # Check plaquettes.
    for i, p in enumerate(builder.plaquettes):
        assert expected[i] == p[:-1]

    # Check the mask array.
    for i, p in enumerate(builder.plaquettes):
        mask = 0
        for k in expected[i]:
            mask = mask + (1 << k)
        assert mask == p[-1]


def test_plaquette_list_3D():
    """ Check if the correct list of plaquettes is found in 3D.
    """
    param = {
        'L' : [2,2,2]
    }
    builder = HamiltonianBuilder(param, states=[])

    # Set up test case.
    expected = [
        [0, 4, 6, 1],
        [1, 8, 13, 2],
        [0, 5, 12, 2],

        [3, 1, 9, 4],
        [4, 11, 16, 5],
        [3, 2, 15, 5],

        [6, 10, 0, 7],
        [7, 2, 19, 8],
        [6, 11, 18, 8],

        [9, 7, 3, 10],
        [10, 5, 22, 11],
        [9, 8, 21, 11],


        [12, 16, 18, 13],
        [13, 20, 1, 14],
        [12, 17, 0, 14],

        [15, 13, 21, 16],
        [16, 23, 4, 17],
        [15, 14, 3, 17],

        [18, 22, 12, 19],
        [19, 14, 7, 20],
        [18, 23, 6, 20],

        [21, 19, 15, 22],
        [22, 17, 10, 23],
        [21, 20, 9, 23],
    ]

    # Chekc the plaquette index.
    for i, p in enumerate(builder.plaquettes):
        assert expected[i] == p[:-1]

    # Check the mask array.
    for i, p in enumerate(builder.plaquettes):
        mask = 0
        for k in expected[i]:
            mask = mask + (1 << k)
        assert mask == p[-1]


def test_u_operator():
    """ Check if the plaquette operator U does what it should do.
    """
    # Testcase 1: should give an overlap.
    # |0000 0001 0010 0000> ---> |0000000010010000>
    builder = HamiltonianBuilder({'L':[2,4]}, states=[])
    state = int('0000000100100000', 2)
    expected_state = int('0000000010010000', 2)
    constructed_state = builder.apply_u(state, builder.plaquettes[2])
    assert expected_state == constructed_state

    # Same test, but slightly different.
    state = int('0000000100100000', 2)
    expected_state = int('0000000010010000', 2)
    constructed_state = HamiltonianBuilder.apply_u(None, state, [4, 7, 8, 5, 432])
    assert expected_state == constructed_state

    # Testcase 2: should anihilate
    # |0010 0000 0010 0000> ---> 0
    state = int('0010 0000 0000 0000'.replace(' ', ''), 2)
    constructed_state = HamiltonianBuilder.apply_u(None, state, [4, 7, 8, 5, 432])
    assert constructed_state == 0

    # ---

    # Testcase 3: 3D builder.
    builder = HamiltonianBuilder({'L':[2,2,2]}, states=[])
    state = _set_bits([20,7])
    expected = _set_bits([19,14])

    p_ind = [19, 14, 7, 20]
    constructed = HamiltonianBuilder.apply_u(None, state, p_ind + [_set_bits(p_ind)])
    assert expected == constructed

    constructed = builder.apply_u(state, builder.plaquettes[19])
    assert expected == constructed



def test_u_dagger_operator():
    """ Check if the inverse plaquette operator U^+ does what it should do.
    """
    builder = HamiltonianBuilder({'L':[2,2,2]}, states=[])
    state = _set_bits([20,7])
    expected = _set_bits([19,14])

    p_ind = [19, 14, 7, 20]
    constructed = HamiltonianBuilder.apply_u_dagger(None, state, p_ind + [_set_bits(p_ind)])
    assert 0 == constructed

    constructed = builder.apply_u_dagger(state, builder.plaquettes[19])
    assert 0 == constructed

    # ---

    state = _set_bits([19,14])
    expected = _set_bits([20,7])

    p_ind = [19, 14, 7, 20]
    constructed = HamiltonianBuilder.apply_u_dagger(None, state, p_ind + [_set_bits(p_ind)])
    assert expected == constructed
