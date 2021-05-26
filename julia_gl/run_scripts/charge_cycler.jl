#===============================================================================

    charge_cycler.jl - LR, May 2021.

===============================================================================#
ENV["JULIA_DEBUG"]="all"
@debug "Debugging activated!"

include("../src/io/io.jl")
include("../src/io/logging.jl")
include("../src/param_checks.jl")
include("snippets/_set_type.jl")


# Read the config file (first argument after the program name).
# param = param_checks!(read_config(ARGS[1]))
raw_config = read_config(length(ARGS)>0 ? ARGS[1] : "config.yml")
make_logger(raw_config)

# Set the size for the Link Type.
const LinkType = set_type(raw_config)

# Include this only after the Type was set.
include("snippets/_multilambda_cycler.jl")

# Get all winding sectors and skip specified ones.
cycle_charges = [
    [[1], [15]],
]

# cycle_charges = [
#     [[1], [2]],
#     [[1], [6]],
#     [[1], [8]],
# ]
# cycle_charges = [[[1], [16]]]

# Go through the loop.
@info "===================== starting ws cycle ========================="
log_error(raw_config, "Diagonalization errors in winding sectors"; init=true)
for static_charges in cycle_charges
    @info "Charge config $static_charges in progress."
    raw_config["static_charges"] = static_charges
    param = copy(raw_config) # Attention: copy is important here!
    # try
        multilambda_cycler(param_checks!(param))
    # catch LoadError
    #     log_error(param, "charge config: $static_charges")
    #     @error "Something went wrong with the diagonalization." static_charges=static_charges
    # end
end
@info "============================== fin =============================="
