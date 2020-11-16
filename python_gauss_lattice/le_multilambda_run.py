""" ----------------------------------------------------------------------------

    le_multilambda_run.py - LR, October 2020

    Scans through parameter space starting from a constructed Hamiltonian for
    low energy states.

---------------------------------------------------------------------------- """
from gl_simulation import LowEnergyGLSimulation
from gauss_lattice import ParallelHamiltonianBuilder
import numpy as np


sim = LowEnergyGLSimulation()


# Produce the Hamiltonian.
grp_name = "ex{:d}".format(sim.param["maximum_excitation_level"])
ham_name = "le_hamiltonian_"+grp_name
ham = sim.read_hamiltonian(ham_name=ham_name)

if not ham:

    # Produce states.
    states, _ = sim.read_le_states(sim.param['maximum_excitation_level'])
    if not len(states):
        base_lattices = {
            (2,2,2) : [3816540, 3872106, 5421780, 5678001, 7542990, 7743645,
                            9033570, 9234225, 11099214, 11355435, 12905109, 12960675],
            (2,2,4) : [64030919769180, 64963162609002, 90962379586260, 95261054903217, 126550380058830, 129916812535965,
                            151558164174690, 154924596651825, 186213921807438, 190512597124395, 216511814101653, 217444056941475],
            (2,2,6) : [1074260571646206819420, 1089901011134353970538,  1526095490192680073940,     1598215294499136381873,
                            2123163061129091160270, 2179642425947400317085, 2542724056922244896610, 2599203421740554053425,
                            3124151188370508831822, 3196270992676965139755, 3632465471735291243157, 3648105911223438394275],
        }
        states = sim.find_le_states(base_lattices[tuple(sim.param['L'])])

    ham = sim.construct_hamiltonian(states, builder_type=ParallelHamiltonianBuilder)
    if sim.store_ham:
        sim.store_hamiltonian(ham, label=ham_name, attrs={'n_fock':len(states)})


# Perform the actual diagonalization.
lambdas = np.linspace(*sim.param['lambdas']) if 'lambdas' in sim.param else [sim.param['lambda']]
sim.run_lambda_loop(lambdas, ham, grp_name=grp_name)
sim.push_notification('Done with multiparameter run.')
