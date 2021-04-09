#===============================================================================

    es_converter.jl - LR, March 2021

    This is a quick check for evaluating expectation values of operators.

===============================================================================#
using DataFrames, CSV
include("../src/typedefs.jl")
include("../src/operators/gl_operators.jl")

# Only for small lattices.
const LinkType = SmallLinkState

function read_states(filename::String)#::Array{WaveFunction,1}
    df = DataFrame(CSV.File(filename, header=0))
end


# Chose operator.
op = gl_operators["2x2x2"]["parity"]
@info "Showing results for operator." name=op.name

# Get Eigenstates.
estates = read_states("/home/lukas/projects/QLMs/FQLMs/parity_troubleshooting/fermiLinks/evecs_lamm100.dat")
print(estates)
#
# for state in states
#     new_state, new_sign = op*state
#     println(Int(state), " ==> ", new_sign==1 ? "+" : "-", Int(new_state))
# end
