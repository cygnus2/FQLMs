import numpy as np
from ed import *

from scipy.sparse import spdiags
from scipy.sparse.linalg import eigsh

L = 4
U = 4
μ = 0

print("L is ", L)
print("U is ", U)
print("μ is ", μ)

H0 = H_free_apbc(L, μ=μ)
H02 = H_free(L, m=1)
HU = H_int(L, U=U)


E0, EU, v0 = energies(L, H02, HU)
print(
    "\nExact ED energies are: \t\t",
    E0,
    "\t",
    EU,
    "sum is ",
    E0 + EU
)