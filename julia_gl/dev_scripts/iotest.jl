include("../src/io/python_import.jl")
include("../src/lattice.jl")

param = read_config(ARGS[1])
(ws, lookup_table, ilookup_table) = read_lookup_tables(param)
