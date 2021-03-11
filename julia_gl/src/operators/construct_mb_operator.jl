#===============================================================================

    construct_mb_operator.jl - LR, January 2021

    This provides functionality to create an operator in a given many-body
    Hilbert space.

===============================================================================#
include("../typedefs.jl")

struct HilbertOperator
    """ Sparse representaiton of an operator that acts on a Hilbert space.
    """
    row::Vector{HilbertIndex}
    col::Vector{HilbertIndex}
    data::Vector{DType}
    n_fock ::HilbertIndex
end

function construct_operator(
        op::Operator,
        lookup_table::LookupDict,
        ilookup_table::InvLookupDict;
        compute_sign::Bool=true
    )::HilbertOperator
    """ Constructs the sparse matrix that represents a given operator in a given
        many-body Hilbert space.

        Note: this is very similar to what is done for the Hamiltonian. Here,
        however, we assume that application of an operator exactly gives one
        state, i.e., the operator consists only of a single term. Multiple terms
        should probably be supproted at some point - right now we don't need this.
    """
    # Loop through all Fock states and create the overlap matrix.
    row = Vector{HilbertIndex}()
    col = Vector{HilbertIndex}()
    data = Vector{CType}()

    # Loop through the list.
    for k=1:length(ilookup_table)

        # This is where we assume only a single term in the operator. We'd need
        # to change some structure to make this more general.
        new_state, sign = apply_operator(op, lookup_table[k]; sign=compute_sign)

        i = get(ilookup_table, new_state, nothing)
        if !isnothing(i)
            push!(row, k)
            push!(col, i)
            push!(data, sign)
        end
    end
    # println(Int.(row))
    return HilbertOperator(row, col, data, length(lookup_table))
end


function apply_operator(ho::HilbertOperator, state::WaveFunction)::WaveFunction
    """ Applies an operator to a wavefunction and returns the new state.
    """
    sparse_ho = sparse(ho.row, ho.col, ho.data, ho.n_fock, ho.n_fock)
    return sparse_ho*state # Potentially large matrix-vector operation.
end

function expectation_value(ho::HilbertOperator, state::WaveFunction)
    """ Computes the expectation value of a given operator with respect to the
        specified state.

        Notes:
         - The dot product already takes care of complex conjugation.
    """
    return dot(state, apply_operator(ho, state))
end
