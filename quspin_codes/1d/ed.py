import sys
import numpy as np
from quspin.operators import hamiltonian
from quspin.basis import spinful_fermion_basis_1d
from scipy.sparse.linalg import expm_multiply
from scipy.optimize import minimize, differential_evolution, newton, basinhopping

np.set_printoptions(precision=5, suppress=True)


def H_free_apbc(L, μ=0):
    basis = spinful_fermion_basis_1d(L)
    hop_right = [[+1, i, i + 1] for i in range(L - 1)]  # OBC
    hop_left = [[-1, i, i + 1] for i in range(L - 1)]  # OBC
    hop_right.append([-1, L - 1, 0])  # APBC
    hop_left.append([+1, L - 1, 0])  # APBC
    pot = [[-μ, i] for i in range(L)]  # -\mu \sum_j n_{j \sigma}I
    # print(hop_right)
    static = [
        ["+-|", hop_left],  # up hops left
        ["|+-", hop_left],  # down hops left
        ["-+|", hop_right],  # up hops right
        ["|-+", hop_right],  # down hops right
        ["n|", pot],  # up on-site potential
        ["|n", pot],  # down on-site potential
    ]
    dynamic = []
    no_checks = dict(check_pcon=False, check_symm=False, check_herm=False)
    H = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, **no_checks)
    return H


def H_free(L, m=1):
    basis = spinful_fermion_basis_1d(L)
    hop_right = [[+1, i, (i + 1) % L] for i in range(L)]  # PBC
    hop_left = [[-1, i, (i + 1) % L] for i in range(L)]  # APBC
    static = [
        ["+-|", hop_left],  # up hops left
        ["|+-", hop_left],  # down hops left
        ["-+|", hop_right],  # up hops right
        ["|-+", hop_right],  # down hops right
    ]
    dynamic = []
    no_checks = dict(check_pcon=False, check_symm=False, check_herm=False)
    H = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, **no_checks)
    return H


def H_int(L, U=4.0):
    basis = spinful_fermion_basis_1d(L)
    interact = [[U, i, i] for i in range(L)]  # U/2 \sm_j n_{j,up} n_{j,down}
    pot = [[-U / 2, i] for i in range(L)]  # -\mu \sum_j n_{j \sigma}I
    static = [
        ["n|n", interact],  # up-down interaction
        ["n|", pot],  # up on-site potention
        ["|n", pot],  # down on-site potent ion
    ]
    dynamic = []
    no_checks = dict(check_pcon=False, check_symm=False, check_herm=False)
    H = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, **no_checks)
    return H


def exact_diag(L, U=4):
    H = H_free(L) + H_int(L, U=U)
    E_GS, V_GS = H.eigsh(k=1, which="SA", maxiter=1e10)  # only GS
    v0 = V_GS[:, 0]
    v0 = v0 / np.linalg.norm(v0)
    return E_GS[0], v0


def energies(L, H0, HU):
    E, v = (H0 + HU).eigsh(k=1, which="SA", maxiter=1e10)  # only GS
    v = v / np.linalg.norm(v)
    E0 = H0.matrix_ele(v, v)
    EU = HU.matrix_ele(v, v)
    return E0, EU, v
