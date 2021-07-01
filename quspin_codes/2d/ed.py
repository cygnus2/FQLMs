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

    J=1.0

    s = np.arange(N_2d)  # sites [0,1,2,....]
    x = s % Lx  # x positions for sites
    y = s // Lx  # y positions for sites
    T_x = (x + 1) % Lx + Lx * y  # translation along x-direction
    P_x = x + Lx * (Ly - y - 1)  # reflection about x-axis
    T_y = x + Lx * ((y + 1) % Ly)  # translation along y-direction
    P_y = (Lx - x - 1) + Lx * y  # reflection about y-axis
    S = -(s + 1)  # fermion spin inversion in the simple case
    P_xy = (Lx - x - 1) + Lx * (Ly - y - 1) # point reflection about origin

    #basis = spinful_fermion_basis_general(
     #   N_2d, pxyblock=(P_xy, 1)
    #)

    basis = spinful_fermion_basis_general(
        N_2d, double_occupancy=False
    )

    #print(basis)

    hop_right = (
        [[+0.5, i, T_x[i]] for i in range(N_2d)]
        + [[+0.5, i, T_y[i]] for i in range(N_2d)]
    )
    hop_left = (
        [[-0.5, i, T_x[i]] for i in range(N_2d)]
        + [[-0.5, i, T_y[i]] for i in range(N_2d)]
    )
    print(hop_right)

    coupling_1 = (
        [[-J/2, i, T_x[i]] for i in range(N_2d)]
        + [[-J/2, i, T_y[i]] for i in range(N_2d)]
        + [[-J/2, T_x[i], i] for i in range(N_2d)]
        + [[-J/2, T_y[i], i] for i in range(N_2d)]
    )

    coupling_2 = (
        [[J/2, i, i, T_x[i]] for i in range(N_2d)]
        + [[J/2, i, i, T_y[i]] for i in range(N_2d)]
        + [[J/2, T_x[i], i, T_x[i]] for i in range(N_2d)]
        + [[J/2, T_y[i], i, T_y[i]] for i in range(N_2d)]
    )

    coupling_3 = (
        [[J/2, i, T_x[i], i] for i in range(N_2d)]
        + [[J/2, i, T_y[i], i] for i in range(N_2d)]
        + [[J/2, i, T_x[i], T_x[i]] for i in range(N_2d)]
        + [[J/2, i, T_y[i], T_y[i]] for i in range(N_2d)]
    )

    coupling_4 = (
        [[-J, i, T_x[i], i, T_x[i]] for i in range(N_2d)]
        + [[-J, i, T_y[i], i, T_y[i]] for i in range(N_2d)]
    )

    coupling_5 = (
        [[-J/2, i, T_x[i], T_x[i], i] for i in range(N_2d)]
        + [[-J/2, i, T_y[i], T_y[i], i] for i in range(N_2d)]
        + [[-J/2, T_x[i], i, i, T_x[i]] for i in range(N_2d)]
        + [[-J/2, T_y[i], i, i, T_y[i]] for i in range(N_2d)]
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

    num = [[1, i] for i in range(N_2d)]

    particle_number = [
        ["n|", num],
        ["|n", num],
    ]

    particle_up = [
        ["n|", num],
    ]

    tJ = [
        ["n|n", coupling_1],
        ["n|nn", coupling_2],
        ["nn|n", coupling_3],
        ["nn|nn", coupling_4],
        ["+-|+-", coupling_5]
    ]

    return basis, static_free, static_int, slater_free, particle_number, particle_up, tJ
