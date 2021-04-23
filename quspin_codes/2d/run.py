import numpy as np
from ed import *
from scipy.sparse import spdiags
from scipy.sparse.linalg import eigsh

Lx = 4
Ly = 2
U = 2
μ = 0

#print("\nL is ", L)
print("U is ", U)
print("μ is ", μ)

basis, free, interact, slater_free, particle_number, particle_up = build_basis(Lx, Ly, U=U, μ=μ)

H0 = hamiltonian(free, [], basis=basis, dtype=np.float64, check_symm=False, check_pcon=False)
HU = hamiltonian(interact, [], basis=basis, dtype=np.float64, check_symm=False, check_pcon=False)
Nop = hamiltonian(particle_number, [], basis=basis, dtype=np.float64, check_symm=False, check_pcon=False)
Nup = hamiltonian(particle_up, [], basis=basis, dtype=np.float64, check_symm=False, check_pcon=False)


E, psi0 = (H0 + HU).eigsh(k=1, which="SA", maxiter=1e6)
print("\nExact E is :\t", E)
print("\n Lowest state is:\t")
#print(psi0[:,0])
print("\n Particle Number is:\t",Nop.expt_value(psi0,time=0))
print("\n Number of up particles is:\t",Nup.expt_value(psi0,time=0))