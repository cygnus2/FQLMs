#===============================================================================

    hamiltonian.jl - LR, January 2021

    Holds stuff for the definition of the Hamiltonian.

===============================================================================#
using SparseArrays, Arpack


struct GaussLatticeHamiltonian <: Hamiltonian
    row::Vector{HilbertIndex}
    col::Vector{HilbertIndex}
    data::Vector{DType}
    n_fock ::HilbertIndex

    # GaussLatticeHamiltonian(row::Vector{HilbertIndex},col::Vector{HilbertIndex},data::Vector{DType}) = begin
    #     if (length(data) != length(row)) || (length(data) != length(row))
    #         error("Faulty data speficifed for sparse Hamiltonian.")
    #     end
    #     new(row,col,data,length(data))
    # end
end

function diagonalize(
    ham::GaussLatticeHamiltonian,
    n_eigenvalues::Integer,
    ev_type::String,
    J::DType,
    lambda::DType,
    gp::String
)
    scaled_data = (gp == "bosons") ? J.*abs.(ham.data) : J.*ham.data
    sparse_ham = sparse(ham.row, ham.col, scaled_data, ham.n_fock, ham.n_fock)

    # If an interaction term exists, we must construct it.
    # if abs(lam):
    #     diag = np.zeros(self.n_fock)
    #     for k in self.row:
    #         diag[int(k)] += lam
    #     self.sparse_rep.setdiag(diag)

    # Perform the usual diagonalization.
    which_map = Dict([
        ("SA", :SR)
    ])

    return eigs(sparse_ham, nev=n_eigenvalues, which=which_map[ev_type])
    # ev, est =  eigs(sparse_ham, nev=n_eigenvalues, which=:SR)
end
