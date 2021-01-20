include("../src/plaquettes.jl")

p = Plaquette(Array{UInt8,1}([1, 2, 3, 5]))
for l in p
    println(l)
end
println(p.mask)
