import matplotlib.pyplot as plt
from lattice_object2D import LatticeObject2D



latt = LatticeObject2D(9248312981713987977913, [6,6])

fig, ax = plt.subplots()
fig.set_size_inches(5,5)
ax.set_aspect(1)
latt.draw(axis=ax, type='spin')
fig.savefig('lattice.png')
