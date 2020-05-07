""" ----------------------------------------------------------------------------

    hamiltonian.py - LR, May 2020

    Hamiltonian object for easier handling.

---------------------------------------------------------------------------- """
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh, eigs


class Hamiltonian(object):
    """ Represents the Hamiltonian.
    """
    def __init__(self, data, row, col, shape):
        """ Takes is a sparse matrix that holds the entries of the hamiltonian.
        """
        self.sparse_rep = csr_matrix(
            (data, (row, col)),
            shape=shape,
            dtype=np.float
        )
        self.diagonalized = False

    def __repr__(self):
        return str(self.sparse_rep.todense())


    def compute_lower_spectrum(self, n_eigenvalues=20, which='LM', dense=False):
        """ Performs the diagonalization with ARPACK and returns the lower part
            (n_eigenvalues) of the spectrum.

            Relies on the method scipy.sparse.linalg.eigsh, which is able to
            treat hermitean matrices (such as the Hamiltonian).
        """
        self.diagonalized = True
        if dense:
            # self.eigenvalues, _ = eigs(self.sparse_rep.todense(), n_eigenvalues)
            self.eigenvalues, _ = eigsh(self.sparse_rep.todense(), n_eigenvalues, which=which)
        else:
            # self.eigenvalues, _ = eigs(self.sparse_rep, n_eigenvalues)
            self.eigenvalues, _ = eigsh(self.sparse_rep, n_eigenvalues, which=which)

        self.eigenvalues = sorted(self.eigenvalues)
        return self.eigenvalues
        # raise NotImplementedError('Not done yet - sorry.')


    def store_results(self, filename='diagonalization_results.dat'):
        """ Dumps the results to file.
        """
        if self.diagonalized:
            with open(filename, 'w') as f:
                for e in self.eigenvalues:
                    f.write('{:.8f},{:.8f}\n'.format(e.real, e.imag))
        else:
            raise ValueError('Hamiltonian cannot be stored yet - it is not diagonalized!')
