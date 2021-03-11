#===============================================================================

    op_lattice_check.jl - LR, March 2021

    Checks the parity for a set of states.

===============================================================================#
include("../src/typedefs.jl")
include("../src/operators/gl_operators.jl")

# Tests only for small lattices.
const LinkType = SmallLinkState

# States to check. Those are the 12 most flippable states for a 2x2x2 lattice.
states = LinkType.([3816540, 3872106, 5421780, 5678001, 7542990, 7743645, 9033570, 9234225, 11099214, 11355435, 12905109, 12960675])

# Sanity: this must map to -9.
# states = LinkType.([9])

# Chose operator.
op = gl_operators["2x2x2"]["parity"]
@info "Showing results for operator." name=op.name
for state in states
    new_state, new_sign = op*state
    println(Int(state), " ==> ", new_sign==1 ? "+" : "-", Int(new_state))
end
