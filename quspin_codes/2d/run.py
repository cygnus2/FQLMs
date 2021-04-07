import numpy as np
from ed import *
from scipy.sparse import spdiags
from scipy.sparse.linalg import eigsh

L = 2
U = 2
μ = 0

print("\nL is ", L)
print("U is ", U)
print("μ is ", μ)

basis, free, interact, slater_free = build_basis(L, L, U=U, μ=μ)

H0 = hamiltonian(free, [], basis=basis, dtype=np.float64, check_symm=False, check_pcon=False)
HU = hamiltonian(interact, [], basis=basis, dtype=np.float64, check_symm=False, check_pcon=False)

E, psi0 = (H0 + HU).eigsh(k=1, which="SA", maxiter=1e6)
print("\nExact E is :\t", E[0])