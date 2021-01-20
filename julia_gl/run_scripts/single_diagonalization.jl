#===============================================================================

    single_diagonalization.jl - LR, January 2021

    Serves to diagonalize a single value. First test with Julia. In essence,
    this is to be replaced with a multilambda run later.

===============================================================================#
using Logging, LoggingExtras
include("../src/io/python_import.jl")
# Read the config file (first argument after the program name).
param = read_config(ARGS[1])

# Set up a logger that persists to file and writes to the console.
logger = TeeLogger(
    ConsoleLogger(),
    SimpleLogger(open(haskey(param, "logfile") ? param["logfile"] : "logfile.log", "w")),
)
global_logger(logger)

# First, read the lookup tables (+ inverse lookup table).
(ws, lookup_table, ilookup_table) = read_lookup_tables(param)

# Then, construct the Hamiltonian.
# TODO.

# Finally, diagonalize.
# TODO.
