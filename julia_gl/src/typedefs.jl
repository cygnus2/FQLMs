#===============================================================================

    typedefs.jl - LR, January 2021

    Holds some definitions of types.

===============================================================================#
# Base type aliases for constant usage.
const DType = Float64
const CType = Complex{DType}
const IType = Int64

const SiteIndex = UInt8
const PlaquetteIndex = UInt8
const HilbertIndex = UInt32
const WaveFunction = Array{CType,1}

# ---
# Representation of the problem.

# Choosing the files actually chooses the implementations - this should only be
# done once, i.e., here. The interface must be the same for all representations.

# LinkState.
include("latt_structure/state.jl")

# Vertex.
# (needs to be after the LinkState definition)
include("latt_structure/vertex.jl")

# Plaquette.
# (needs to be after the LinkState definition)
include("latt_structure/plaquette.jl")

# LinkLattice.
# (needs to be after the definitions of LinkState, Vertex and Plaquette).
abstract type Lattice end
include("latt_structure/lattice.jl")

# ---
# Operator stuff.
# abstract type Operator end
include("operators/operators.jl")

# ---
# Diagonalization / construction purposes.
abstract type Hamiltonian end
include("gl_hamiltonian.jl")

# Lookup dictionaries.
const LookupDict = Dict{HilbertIndex,LinkState}
const InvLookupDict = Dict{LinkState,HilbertIndex}
