Python implementation for the Gauss lattice
================================================================================
These are a few scripts that deal with spins/fermions as the quantum links on a
pyrochlore lattice while obeying some sort of Gauss' law. The implementation
should work fine for both 2D and 3D systems.

Generally there a two parts:
 - finding the Gauss law states
 - constructing and diagonalizing the Hamiltonian

Prerequisites:
 - Python 3 (numpy, scipy, h5py, pyyaml, tdqm)

All script, which are described below, can be called with `python <scriptname> -h`, which displays a brief help message that should clarify the details.

 ## Finding the states
The code to find states on a d-dimensional lattice that obey Gauss' Law (based on recursion) is contained in `gauss_lattice/gauss_lattice.py`. The script to run

Run with
```
  python state_finder.py -i <config_file>
```
where the `<config_file>` should be in YAML format. An example would look like this:
```
  L : [2,2,2]  
  working_directory: "./"
```
There can be arbitrarily many entries in the file (ordering is not important) but `L` and `working_directory` need to be specified in order for the script `state_finder.py` to work.

The script will produce (if it not exists) a directory at the location specified by `working_directory` where you can find files containing the number of states per winding sector (`'winding_sectors_<LxLyLz>.dat'`) as well as all the states themselves (`'winding_states_<LxLyLz>.hdf5'`), represented by a sorted list of integers.


## Construction and diagonalization of the Hamiltonian (single value of lambda)
To perform the actual diagonalization of the Hamiltonian run
```
  python diagonalizer.py
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
```

The script produces a bunch of files:
 - `hamiltonian_<LxLyLz>.npz` (holds the Hamiltonian for later use - saves some effort)
 - `spectrum_<gauge_particles>_<LxLyLz>_lam<lambda>.dat.npy` (the specified `n_eigenvalues` of the spectrum, sorted)
 - `eigenstates_<gauge_particles>_<LxLyLz>_lam<lambda>.dat.npy` (holds the eigenstates, only if `compute_eigenstates` is provided and set to `True`)

This produces a file that contains the lower spectrum of the Hamiltonian. In the script, several parameters may be adjusted - the only requirement is, that the `state_finder.py` script was executed for the lattice configuration, otherwise no states can be read from file.

### Sequential diagonalization
This is a special case


## Construction and diagonalization of the Hamiltonian (multiple values of lambda)
