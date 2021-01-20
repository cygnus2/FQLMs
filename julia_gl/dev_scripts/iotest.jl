include("../src/io/python_import.jl")
include("../src/lattice.jl")

param = read_config(ARGS[1])
latt = LinkLattice(param["L"])
ws = haskey(param, "winding_sector") ? _winding_tag(param["winding_sector"], latt=latt) : "all-ws"

s = read_states(param, ws)
println(length(s))
println(typeof(s))
