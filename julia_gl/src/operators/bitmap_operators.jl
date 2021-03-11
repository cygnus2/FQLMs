#===============================================================================

    operators.jl - LR, January 2021

    Represents an operator on a *single lattice configuration*. This essentially
    holds the map where indicies are mapped.

===============================================================================#
using Combinatorics
using LinearAlgebra

abstract type Operator end


struct BitmapOperator <: Operator
    name::String # Just to be explicit.
    index_map::Vector{LinkIndex} # List where indicies are mapped.
end

function Base.show(io::IO, op::BitmapOperator)
    print(io, op.name*" operator (bitmap, length=$(length(op.index_map)))")
end

function LinearAlgebra.dot(op::BitmapOperator, state::LinkState)::Array{SiteIndex}
    """ Returns the indicies of the new particle positions, usorted - not in
        normal order. Those can be used to infer the sign and the new state.

        Example: If the action of the operator is OP[00101] = [11000] then this
        function returns [4, 5] since those are the new positions (in the appropriate
        order as given in the index map of the operator).

        Current implementation: works with repeaded calls to getindex(state),
        which could likely be optimized by updating a state and doing incremental
        single bit shifts.
    """
    # Note: simply reversing this does not just alter a global sign of the
    # permutation. For instance:
    #   sign([1,2,3,4]) = sign([4,3,2,1])
    #   but sign([1,2,3]) = -sign([3,2,1])
    # Interestingly though, reversing does not change the expectation values.
    shuffled = [op.index_map[k] for k=1:length(op.index_map) if state[k]==1]
    return shuffled
end

function apply_operator(op::BitmapOperator, state::LinkState; sign::Bool=true)::Tuple{LinkState,IType}
    """ Applies a bitmap operator and also computes the sign.
    """
    # First, get the transformed (shuffled) list of indicies and produce the state.
    shuffled = dot(op, state)
    new_state = LinkType(0)
    for s in shuffled
        # Shift the indicies, note, those are 1 based.
        new_state += (LinkType(1)<<(s-1))
    end

    if !sign
        return new_state, 1
    end

    # Compute sign by parity of permutation of the sorted state.
    # Note: sorting needs to be ascending, since this is they way the dot product
    # returns the shuffled indicies.
    sign = iseven(parity(sortperm(shuffled))) ? 1 : -1
    return new_state, sign
end

# Provide shorthand notation.
Base.:*(op::BitmapOperator, state::LinkState)::Tuple{LinkState,IType} = apply_operator(op,state;sign=true)
