include("../io/python_import.jl")

param = read_config(ARGS[1])
s = read_states(param, "wx_2-wy_2-wz_2")
println(length(s))
println(typeof(s))
