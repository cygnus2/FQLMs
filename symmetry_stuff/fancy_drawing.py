from lattice_object import LatticeObject


latt = LatticeObject(34232342222342, [2,2,2])
latt.apply_charge_conjugation()
latt.apply_charge_conjugation()
latt.draw()
latt.fig.savefig('fancy_grid.pdf', density=600)
