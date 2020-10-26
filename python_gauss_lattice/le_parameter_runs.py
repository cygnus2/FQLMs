""" ----------------------------------------------------------------------------

    parameter_runs.py - LR, May 2020

    Scans through parameter space starting from a constructed Hamiltonian.

---------------------------------------------------------------------------- """
from gauss_lattice import GaussLatticeHamiltonian, LowEnergyHamiltonianBuilder
from gauss_lattice.aux import size_tag, timeit, read_all_states, file_tag, load_config, read_winding_sector
import numpy as np
import h5py as hdf
import argparse, logging, os

parser = argparse.ArgumentParser(description="Python gauss lattice diagonalizer.")
parser.add_argument('-i', metavar='', type=str, default=None, help='YAML style input file.')
parser.add_argument('-readonly', action='store_true', help='If this is set, the script only runs if the necessary files are found in the working directory (good for cluster usage if large files would be produced).')
parser.add_argument('-notify', action='store_true', help='Toggles whether a notification should be sent when done. Should only be used by Lukas (sorry).')
args = parser.parse_args()

def read_le_states(data_file, max_level, combine=True):
    """ Reads all states up to a certain level.
    """
    states = {}
    with hdf.File(data_file, 'r') as f:
        levels = []
        for g in f:
            l = int(g.split("_")[-1])
            levels.append(l)
            if l <= max_level:
                states[l] = f[g][...]

        if sorted(levels)[-1] < max_level:
            print(f"Warning: maximal level not reachable from list! Taking maximum of {l}.")
            new_max = l
        else:
            new_max = max_level

    if combine:
        return np.concatenate([states[k] for k in sorted(states.keys()) if len(states[k])]), new_max
    return states, new_max


# This sets the parameters fof the calculation (everything else is fixed).
param = load_config(args.i)
hamiltonian_file = param['working_directory'] + '/LE_hamiltonian_' + size_tag(param['L']) + '.hdf5'

# Set up a logger with a handler for the terminal output.
logger = logging.getLogger('parameter run logger')
logger.addHandler(logging.StreamHandler())
logger.addHandler(logging.FileHandler(param['working_directory'] + '/' + param['logfile'], mode='w'))
logger.setLevel(logging.DEBUG)


@timeit(logger=logger)
def hamiltonian_diagonalization(ham, **kwargs):
    """ This method just exists to be able to time the routine efficiently.
    """
    return ham.diagonalize(**kwargs)

@timeit(logger=logger)
def hamiltonian_construction(builder, *args, **kwargs):
    """ This method just exists to be able to time the routine efficiently.
    """
    return builder.construct(*args, **kwargs)

base_lattices = {
    (2,2,2) : [3816540, 3872106, 5421780, 5678001, 7542990, 7743645,
                9033570, 9234225, 11099214, 11355435, 12905109, 12960675],
    (2,2,4) : [64030919769180, 64963162609002, 90962379586260, 95261054903217, 126550380058830, 129916812535965,
                151558164174690, 154924596651825, 186213921807438, 190512597124395, 216511814101653, 217444056941475],
    (2,2,6) : [1074260571646206819420, 1089901011134353970538,  1526095490192680073940, 1598215294499136381873,
                2123163061129091160270, 2179642425947400317085, 2542724056922244896610, 2599203421740554053425,
                3124151188370508831822, 3196270992676965139755, 3632465471735291243157, 3648105911223438394275],
}


# Start by constructing the builder.
le_builder = LowEnergyHamiltonianBuilder(param, logger=logger)

# Check if states exist, if not, construct them.
state_file = param.get("state_file")
if not state_file:
    state_file = param['working_directory']+file_tag(param['L'], filetype='hdf5').replace("winding_", "le_")
if not os.path.isfile(state_file):
    logger.info('Could not find stored states, constructing low-energy Fock state list from scratch.')
    states = le_builder.find_all_states(
        base_lattices[tuple(param['L'])],
        n_threads=param.get('n_threads', 1),
        max_level=param.get('maximum_excitation_level', 10000)
    )
    logger.info('Constructed states.')
else:
    logger.info(f'Reading Fock states from {state_file}')
    states, new_max = read_le_states(state_file, param["maximum_excitation_level"])
    le_builder.read_le_states(states, from_file=True)

# ---
# Retrieves the Hamiltonian and if necessary, constructs it.
ham_name = "le_hamiltonian_ex{:d}".format(param["maximum_excitation_level"])
try:
    with hdf.File(hamiltonian_file, 'r') as f:
        mat = f[ham_name][...]
        ham = GaussLatticeHamiltonian(mat[2,:], mat[1,:], mat[0,:], len(states))
    logger.info('Read Hamiltonian from provided file.')

except (KeyError, OSError):
    # Set up the builder object & construct the Hamiltonian.
    ham = hamiltonian_construction(le_builder, param.get('n_threads', 1))

    # If specified, store the Hamiltonian for later use.
    if param['store_hamiltonian']:
        with hdf.File(hamiltonian_file, 'a') as f:
            if ham_name in f:
                del f[ham_name]
            ds = f.create_dataset(ham_name, data= np.array([ham.col, ham.row, ham.data]))
            ds.attrs['n_fock'] = ham.n_fock

# ----
# Loop through spectra.
spectrum_file = (
    param['working_directory'] +
    '/LE_multi_spectrum_' +
    param['gauge_particles'] + '_' +
    size_tag(param['L']) +
    '.hdf5'
)

lambdas = np.linspace(*param['lambdas'])
for i, l in enumerate(lambdas):
    logger.info('[{:d} / {:d}] diagonalizing Hamiltonian for lambda={:.4f}'.format(i+1, len(lambdas), l))

    # Diagonalization.
    spectrum = hamiltonian_diagonalization(ham,
        full_diag = param.get('full_diag'),
        J = param['J'],
        lam = l,
        gauge_particles = param['gauge_particles'],
        n_eigenvalues = max(1, min(param['n_eigenvalues'], ham.n_fock//2)),
        which = param['ev_type']
    )

    # Output.
    with hdf.File(spectrum_file, 'a') as f:
        f.attrs['lambdas'] = lambdas
        grp_name = "ex_{:d}".format(param["maximum_excitation_level"])
        if grp_name in f:
            grp = f[grp_name]
        else:
            grp = f.create_group(grp_name)
            grp.attrs['maximum_excitation_level'] = param["maximum_excitation_level"]

        # Overwrite dataset.
        ds_name = 'lam_{:6f}'.format(l)
        if ds_name in grp:
            del grp[ds_name]
        ds = grp.create_dataset(ds_name, data=spectrum)
        ds.attrs['lambda'] = l


if args.notify:
    from gauss_lattice.lr_notify import push_message
    import os, subprocess
    host = subprocess.check_output(['hostname']).strip().decode('UTF-8'),
    push_message(f'Multi-parameter run  @{host[0]} is done!')
