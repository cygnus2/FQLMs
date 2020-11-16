from lattice_object import LatticeObject


latt = LatticeObject(34232342222342, [2,2,2])
latt.draw()
latt.fig.savefig('fancy_grid.pdf', density=600)
