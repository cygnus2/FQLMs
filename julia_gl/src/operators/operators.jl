#===============================================================================

    operators.jl - LR, January 2021

    Represents an operator on a *single lattice configuration*. This essentially
    holds the map where indicies are mapped.

===============================================================================#
struct Operator
    name::String # Just to be explicit.
    index_map::Vector{LinkIndex} # List where indicies are mapped.
end

apply_operator(op::Operator, state::LinkState)::LinkState = state[op.index_map]


function Base.show(io::IO, op::Operator)
    print(io, op.name*" operator ($(length(op.index_map)))")
end
