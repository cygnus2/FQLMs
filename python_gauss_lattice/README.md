Python implementation for the Gauss lattice
================================================================================
These are a few scripts that deal with spins/fermions as the quantum links on a
pyrochlore lattice while obeying some sort of Gauss' law. The implementation
should work fine for both 2D and 3D systems.

Generally there a few use cases:
 - finding the Gauss law states (this should be done first)
 - constructing and diagonalizing the Hamiltonian (for a single datapoint)
 - batch diagonalization (i.e., multiple values of the coupling lambda)

Prerequisites:
 - Python 3 (numpy, scipy, h5py, pyyaml, tdqm)

All scripts, which are described below, can be called with `python <scriptname> -h`, which displays a brief help message that should clarify the details.

 ## Finding the states
The code to find states on a d-dimensional lattice that obey Gauss' Law (based on recursion) is contained in `gauss_lattice/gauss_lattice.py`. Run the script with
```
  python state_finder.py -i <config_file>
```
where the `<config_file>` should be in YAML format. A minimal example would look like this:
```
  L : [2,2,2]  
  working_directory: "./"
```
There can be arbitrarily many entries in the file (ordering is not important) but `L` and `working_directory` need to be specified in order for the script `state_finder.py` to work.

The script will produce (if it not exists) a directory at the location specified by `working_directory` where you can find files containing the number of states per winding sector (`'winding_sectors_<LxLyLz>.dat'`) as well as all the states themselves (`'winding_states_<LxLyLz>.hdf5'`), represented by a sorted list of integers.


## Construction and diagonalization of the Hamiltonian (single value of lambda)
To perform the actual diagonalization of the Hamiltonian run
```
  python diagonalizer.py -i <config_file>
```
There are several (self explanatory) parameters that are read from the `<config_file>`
```
  L : [2,2,2]  
  working_directory: "./"

  gauge_particles : 'bosons' # Either bosons or fermions.
  full_diag : False # If True, the full matrix is diagonalized - extremely costly.
  ev_type : 'BE' # ARPACK style specification of the types of eigenvalues.
  n_eigenvalues : 50 # Number of eigenvalues to compute.
  lambda: -5 # lambda Coupling.
  J : 0.0 # J Coupling.
  compute_eigenstates: True
  # winding_sector : [0,0] # If this is specified, only the appropriate winding sector will be diagonalized.
  compute_eigenstates: True # If true, the eigenstates will be exported (should be True).
  store_hamiltonian: True # If true, the Hamiltonian will be stored (should be True).
```
The script produces a bunch of files:
 - `hamiltonian_<LxLyLz>.npz` (holds the Hamiltonian for later use - saves some effort)
 - `spectrum_<gauge_particles>_<LxLyLz>_lam<lambda>.dat.npy` (the specified `n_eigenvalues` of the spectrum, sorted)
 - `eigenstates_<gauge_particles>_<LxLyLz>_lam<lambda>.dat.npy` (holds the eigenstates, only if `compute_eigenstates` is provided and set to `True`)

The requirement for this script to run is that the `state_finder.py` script was executed for the apropriate lattice configuration, otherwise no states can be read from file - this will result in an I/O error.

### Sequential diagonalization
This is a special case for the diagonalization. The use-case is exactly as above, however, the Hamiltonian will be diagonalized in every winding-sector separately (useful for larger runs). Simply run the script
```
  python sequential_diagonalizer.py -i <config_file>
```
The script works without running anything in advance (except `state_finder.py`) - if the files are not found, they will be generated. In addition, following files will be created:
 - `SEQUENTIAL_spectrum_<gauge_particles>_<LxLyLz>.dat` (Holds the sorted spectrum, combined from all winding sectors)
 - `SEQUENTIAL_spectrum_<gauge_particles>_<LxLyLz>.hdf5` (Holds all spectra for the winding sectors separately in a HDF5 file)
 - `SEQUENTIAL_hamiltonian_<gauge_particles>.hdf5` (Holds the Hamiltonian for each winding sector separately - useful for a multi-parameter run [not implemented right now])
 - `<logfile>` (holds output, filename as specified under `logfile` in the config file - mainly for debugging purposes)


## Construction and diagonalization of the Hamiltonian (multiple values of lambda)
Finally, there's a script that allows the batch diagonalization of several lambda values. Run
```
  python parameter_runs.py -i <config_file>
```
The requirements for the config file are as above, with the exception of the new parameter
```
  lambdas : [-5, 1, 201] # Range that specifies the lambda values to scan (min, max, n_steps)
```
which specifies the range of lambda values.

The script works without running anything in advance (except `state_finder.py`)  - if the files are not found, they will be generated. The output is stored in the file `multi_spectrum_<gauge_particles>_<LxLyLz>.hdf5`, which holds a spectrum for every value of lambda.
