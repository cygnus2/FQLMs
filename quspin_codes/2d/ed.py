import sys
import numpy as np
import math
from quspin.operators import hamiltonian
from quspin.basis import spinful_fermion_basis_general
from scipy.sparse.linalg import expm_multiply
from scipy.optimize import minimize, differential_evolution, newton, basinhopping

np.set_printoptions(precision=5, suppress=True)


def build_basis(Lx, Ly, U=4.0, μ=0):
    N_2d = Lx * Ly

    Nup = N_2d
    Ndown = N_2d

    s = np.arange(N_2d)  # sites [0,1,2,....]
    x = s % Lx  # x positions for sites
    y = s // Lx  # y positions for sites
    T_x = (x + 1) % Lx + Lx * y  # translation along x-direction
    P_x = x + Lx * (Ly - y - 1)  # reflection about x-axis
    T_y = x + Lx * ((y + 1) % Ly)  # translation along y-direction
    P_y = (Lx - x - 1) + Lx * y  # reflection about y-axis
    S = -(s + 1)  # fermion spin inversion in the simple case

    basis = spinful_fermion_basis_general(
        N_2d,
    )

    hop_right = (
        [[+0.5, i, T_x[i]] for i in range(N_2d)]
        + [[+0.5, i, T_y[i]] for i in range(N_2d)]
    )
    hop_left = (
        [[-0.5, i, T_x[i]] for i in range(N_2d)]
        + [[-0.5, i, T_y[i]] for i in range(N_2d)]
    )

    pot = [[-U / 2, i] for i in range(N_2d)]  # -\mu \sum_j n_{j \sigma}I

    static_free = [
        ["+-|", hop_left],  # up hops left
        ["|+-", hop_left],  # down hops left
        ["-+|", hop_right],  # up hops right
        ["|-+", hop_right],  # down hops right
    ]

    # Interaction
    interact = [[U, i, i] for i in range(N_2d)]  # U/2 \sm_j n_{j,up} n_{j,down}

    mass1 = [[(-1) ** (i % 2 + math.floor(i / 2)) * 1e-5, i] for i in range(N_2d)]
    if U <= 0:
        mass2 = [[(-1) ** (i % 2 + math.floor(i / 2)) * 1e-5, i] for i in range(N_2d)]
    else:
        mass2 = [[-1*(-1) ** (i % 2 + math.floor(i / 2)) * 1e-5, i] for i in range(N_2d)]

    static_int = [
        ["n|n", interact],  # up-down interaction
        ["n|", pot],  # up on-site potention
        ["|n", pot],  # down on-site potent ion
    ]

    # Proper slater state
    slater_free = [
        ["+-|", hop_left],  # up hops left
        ["|+-", hop_left],  # down hops left
        ["-+|", hop_right],  # up hops right
        ["|-+", hop_right],  # down hops right
        ["n|", mass1],
        ["|n", mass2],
    ]
    return basis, static_free, static_int, slater_free
