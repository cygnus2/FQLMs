Python implementation for the Gauss lattice
================================================================================
These are a few scripts that deal with spins/fermions as the quantum links on a
pyrochlore lattice while obeying some sort of Gauss' law. The implementation
should work fine for both 2D and 3D systems.

There are two parts:
 - finding the Gauss law states
 - constructing and diagonalizing the Hamiltonian

Prerequisites:
 - Python 3

 ## Finding the states
Code to find states on a d-dimensional lattice that obey some sort of Gauss Law.

Run with
```
  python state_finder.py
```
In the script you can easily change the parameters and the I/O behavior via suitable parameters (to be defined at a later point here).


The script will produce (if it not exists) a directory `output` where you can find several files containing the number of states per winding sector as well (if specified) all the states themselves (in an HDF5 or plain text file, represented by integers).


## Construction and diagonalization of the Hamiltonian
This is done in a separate script (partly to save time, partly for historic reasons).

Run with
```
  python diagonalizer.py
```
This produces a file that contains the lower spectrum of the Hamiltonian. In the script, several parameters may be adjusted - the only requirement is, that the `state_finder.py` script was executed for the lattice configuration, otherwise no states can be read from file.
