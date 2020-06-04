""" ----------------------------------------------------------------------------

    hamiltonian.py - LR, May 2020

    Hamiltonian object for easier handling of data pipeline.

---------------------------------------------------------------------------- """
import numpy as np
from scipy.sparse import coo_matrix, save_npz, load_npz
from scipy.sparse.linalg import eigsh, eigs
from scipy.linalg import eigvals, eig
from .aux import write_simple_spectrum
from copy import copy


class Hamiltonian(object):
    """ Represents a relatively general structure for a Hamiltonian (which
        ultimately is only a sparse matrix) and handling some convenient I/O
        business.
    """
    def __init__(self, data, row, col, n_fock):
        """ Takes is a sparse matrix that holds the entries of the hamiltonian.
        """
        self.n_fock = n_fock

        # The storage should be separate, in order to keep the options flexible
        # regarding a change of parameters.
        self.row = np.array(row)
        self.col = np.array(col)
        self.data = np.array(data)

        self.diagonalized = False
        self.sparsified = False

    @classmethod
    def from_scipy_dump(cls, input_file):
        """ Alternate setup with data from file.
        """
        sp = load_npz(input_file)
        return cls(sp.data, sp.row, sp.col, sp.shape[0])

    @classmethod
    def from_hdf5_dump(cls, input_file, ws):
        """ Alternate setup with data from file.
        """
        raise NotImplementedError('HDF5 loading not implemented yet, sorry!')


    def __repr__(self):
        return str(self.sparse_rep.todense())


    def _sparsify(self, entries=None):
        """ Constructs the sparse matrix from the pieces we stored.
        """
        if entries is None:
            self.sparse_rep = coo_matrix(
                (self.data, (self.row, self.col)),
                shape=(self.n_fock, self.n_fock),
                dtype=np.float
            )
        else:
            self.sparse_rep = coo_matrix(
                entries,
                shape=(self.n_fock, self.n_fock),
                dtype=np.float
            )
        self.sparsified = True


    def store_hamiltonian(self, filename='hamiltonian_sparse.npz'):
        """ Stores the Hamiltonian matrix in sparse format.
        """
        if not self.sparsified:
            self._sparsify()
        save_npz(filename, self.sparse_rep)


    def diagonalize(self, n_eigenvalues=50, which='BE', full_diag=False, compute_eigenstates=False):
        """ The diagonalization routine which is called from the outside. If
            full_diag is set to True, a full eigensolver (not ARPACK) will be
            used which could lead to a dramatic loss of performance.
        """
        if not self.sparsified:
            self._sparsify()

        # Perform the actual diagonalization.
        if full_diag:
            self._full_diagonalization(compute_eigenstates=compute_eigenstates)
        else:
            self._compute_lower_spectrum(n_eigenvalues, which)
        self.diagonalized = True

        # Sorting.
        m = np.argsort(self.eigenvalues)
        if compute_eigenstates:
            return self.eigenvalues[m], self.eigenstates[:,m]
        return self.eigenvalues[m]


    def _compute_lower_spectrum(self, n_eigenvalues, which):
        """ Performs the diagonalization with ARPACK and returns the lower part
            (n_eigenvalues) of the spectrum.

            Relies on the method scipy.sparse.linalg.eigsh, which is able to
            treat hermitean matrices (such as the Hamiltonian).
        """
        # Warning: the diagonal of the Hamiltonian is set to 0 - this is not the
        # most general case.
        if self.n_fock == 1:
            self.eigenvalues, self.eigenstates = [0], [0]
        else:
            self.eigenvalues, self.eigenstates = eigsh(self.sparse_rep, n_eigenvalues, which=which)


    def _full_diagonalization(self, compute_eigenstates=False):
        """ Full diagonalization of the Hamiltonian.

            Clearly, this should/only be done for small systems - for larger ones
            this routine will explode in runtime.
        """
        if compute_eigenstates:
            self.eigenvalues, self.eigenstates = eig(self.sparse_rep.todense())
        else:
            self.eigenvalues = eigvals(self.sparse_rep.todense())


    def store_results(self, filename='diagonalization_results.dat', store_eigenvalues=False):
        """ Dumps the results for the spectrum to file.
        """
        if self.diagonalized:
            write_simple_spectrum(self.eigenvalues, filename)
            if store_eigenvalues:
                np.save(filename.replace('spectrum', 'eigenstates'), self.eigenstates)
        else:
            raise ValueError('Hamiltonian cannot be stored yet - it is not diagonalized!')



class GaussLatticeHamiltonian(Hamiltonian):
    """ Represents a specific implementation of a Hamiltonian for the Gauss
        lattice which reads

            H = -J sum_{p} (U_p + U+_p) + lam sum_{p} (U_p + U+_p)**2

        where the sum is over all plaquettes. The first term is a plaquette-
        flipping term, a correlated hop of two fermions. The second part is


        The concretization allows us to provide parameters
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


    def diagonalize(self, J=1, lam=0, gauge_particles='fermions', **kwargs):
        """ Performs the diagonalization. If no keyword arguments are provided,
            the fermionic Hamiltonian with J=1 and lambda = 0 will be computed.

            Alternatively, parameters may be provided to obtain the spectrum at
            the specific parameter choice.
        """
        if gauge_particles == 'bosons':
            self._sparsify(entries=(J*np.abs(self.data), (self.row, self.col)))
        else:
            self._sparsify(entries=(J*self.data, (self.row, self.col)))

        # If an interaction term exists, we must construct it.
        if abs(lam):
            diag = np.zeros(self.n_fock)
            for k in self.row:
                diag[k] += lam
            self.sparse_rep.setdiag(diag)


        # Perform the usual diagonalization.
        return super().diagonalize(**kwargs)
