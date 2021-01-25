#===============================================================================

    hamiltonian.jl - LR, January 2021

    Holds stuff for the definition of the Hamiltonian.

===============================================================================#
using SparseArrays, LinearAlgebra, Arpack
using PyCall
@pyimport scipy.sparse.linalg as pylinalg
@pyimport scipy.sparse as pysparse


struct GaussLatticeHamiltonian <: Hamiltonian
    row::Vector{HilbertIndex}
    col::Vector{HilbertIndex}
    data::Vector{DType}
    n_fock ::HilbertIndex
end

function diagonalize(ham::GaussLatticeHamiltonian, lambda::DType, param::Dict{Any,Any})
    scaled_data = (param["gauge_particles"] == "bosons") ? param["J"].*abs.(ham.data) : param["J"].*ham.data
    sparse_ham = sparse(ham.row, ham.col, scaled_data, ham.n_fock, ham.n_fock)

    # If an interaction term exists, we must construct it.
    if abs(lambda) > 0
        diag = fill(DType(0.0), ham.n_fock)
        for k in ham.row
            diag[HilbertIndex(k)] += lambda
        end
        sparse_ham[diagind(sparse_ham)] = diag
    end

    # Perform the usual diagonalization.
    which_map = Dict([
        ("SA", :SR)
    ])
    return eigs(
        sparse_ham,
        nev=param["n_eigenvalues"],
        which=which_map[param["ev_type"]],
        tol=get(param, "tol", 1e-15),
        maxiter=get(param, "max_iter", 10*ham.n_fock)
    )
end
# Increasing maxiter actually helped in some cases - might also be the case for lambda = -3.0?
# Julia has a matrix-independent iteration number of 300 per default, scipy has it at 10*N where N is the dimension of the Hilbert space. This probably explains why there's an issue.


# function diagonalize(
#     ham::GaussLatticeHamiltonian,
#     n_eigenvalues::Integer,
#     ev_type::String,
#     J::DType,
#     lambda::DType,
#     gp::String
# )
# """ A version of the diagonalization with calls to scipy, which appears to be
#     working for some reason.
# """
#     scaled_data = (gp == "bosons") ? J.*abs.(ham.data) : J.*ham.data
#     sparse_ham = pysparse.csc_matrix((scaled_data, (ham.row .- 1, ham.col .- 1)), shape=(ham.n_fock,ham.n_fock))
#
#     # If an interaction term exists, we must construct it.
#     if abs(lambda) > 0
#         diag = fill(DType(0.0), ham.n_fock)
#         for k in ham.row
#             diag[HilbertIndex(k)] += lambda
#         end
#         sparse_ham.setdiag(diag)
#     end
#
#     # Perform the usual diagonalization.
#     which_map = Dict([
#         ("SA", :SR)
#     ])
#     # return eigs(sparse_ham, nev=n_eigenvalues, which=which_map[ev_type])
#     return pylinalg.eigsh(sparse_ham, k=n_eigenvalues, which="SA")
#     # ev, est =  eigs(sparse_ham, nev=n_eigenvalues, which=:SR)
# end
