#===============================================================================

    single_diagonalization.jl - LR, January 2021

    Serves to diagonalize a single value. First test with Julia. In essence,
    this is to be replaced with a multilambda run later.

===============================================================================#
ENV["JULIA_DEBUG"]="all"
@debug "Debugging activated!"

include("../src/io/io.jl")
include("../src/io/logging.jl")
include("../src/param_checks.jl")
include("snippets/_set_type.jl")


# Read the config file (first argument after the program name).
# param = param_checks!(read_config(ARGS[1]))
param = param_checks!(read_config(length(ARGS)>0 ? ARGS[1] : "config.yml"))
make_logger(param)

# Set the size for the Link Type.
const LinkType = set_type(param)

# Include this only after the Type was set.
include("snippets/_multilambda_cycler.jl")
multilambda_cycler(param)
