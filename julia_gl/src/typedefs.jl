#===============================================================================

    typedefs.jl - LR, January 2021

    Holds some definitions of types.

===============================================================================#
# Types for state handling. The type of representation is chosen by including the
# corresponding file.
include("state_reps/rep_integer.jl")
const Vertex = Vector{UInt8}
const LookupDict = Dict{UInt32,LinkState}
const InvLookupDict = Dict{LinkState,UInt32}

# Types for storage of numerical coefficients.
const DType = Float64
