import h5py as hdf
import numpy as np

# Read GS.
# datafile = '/home/lukas/projects/qlm_diag/_testing/results_fermions_all-ws_2x2x2.hdf5'
# datafile = '/home/lukas/_TEMP_fqlm/results_fermions_all-ws_2x2x2.hdf5'
datafile = '/home/lukas/_TEMP_fqlm//results_fermions_wx_2-wy_2-wz_2_2x2x2.hdf5'
ds = 'eigenstates_lam_-3.000000'
with hdf.File(datafile, 'r') as f:
    data = f[ds][...][0]

# Read the winding states.
def read_states(filename, ws=None):
    states = []
    with hdf.File(filename, 'r') as f:
        if ws is None:
            for g in f:
                states += list(f[g][...])
        else:
            states = list(f[ws][...])
    return np.sort(states)

states = read_states('/home/lukas/projects/QLMs/FQLMs/python_data/local_state_storage/winding_states_2x2x2.hdf5', ws='wx_2-wy_2-wz_2')
contributions = [3816540, 3872106, 5421780, 5678001, 7542990, 7743645, 9033570, 9234225, 11099214, 11355435, 12905109, 12960675]

indicies = [np.argwhere(states==c)[0][0] for c in contributions]
print(indicies)
np.savetxt('gs_fermions_2x2x2.txt', data[indicies])


import sys
sys.path.append("../symmetry_stuff/")
from lattice_object import LatticeObject

for i in range(6):
    latt = LatticeObject(contributions[i], [2,2,2])
    latt.apply_charge_conjugation()
    print(latt.to_int(), latt.to_int() == contributions[-i-1])
    print("--")


state = data
print("CC: ", np.dot(state, state[::-1]))
