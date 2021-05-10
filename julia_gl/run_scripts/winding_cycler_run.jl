#===============================================================================

    winding_cycler_run.jl - LR, January 2021.

===============================================================================#
ENV["JULIA_DEBUG"]="all"
@debug "Debugging activated!"

include("../src/io/io.jl")
include("../src/io/logging.jl")
include("../src/param_checks.jl")
include("snippets/_set_type.jl")


function get_all_winding_sectors(L::Array{Int,1})::Set{Tuple}
    if length(L) == 2
        error("not implemented for 2D")
    elseif length(L) == 3
        nx, ny, nz = Int(L[2]*L[3]/2), Int(L[1]*L[3]/2), Int(L[2]*L[1]/2)

        all_ws = Set([])
        for wx=-nx:nx
            for wy=-ny:ny
                for wz=-nz:nz
                    push!(all_ws, (wx,wy,wz))
                end
            end
        end
        return all_ws
    end
    error("Some wrong lattice dimension were specified.")
end


# Read the config file (first argument after the program name).
# param = param_checks!(read_config(ARGS[1]))
raw_config = read_config(length(ARGS)>0 ? ARGS[1] : "config.yml")
make_logger(raw_config)

# Set the size for the Link Type.
const LinkType = set_type(raw_config)

# Include this only after the Type was set.
include("snippets/_multilambda_cycler.jl")

# Get all winding sectors and skip specified ones.
skip_ws = Set([(0,0,0),(1,0,0),(0,1,0),(0,0,1),(-1,0,0),(0,-1,0),(0,0,-1)])
all_ws = get_all_winding_sectors(raw_config["L"])
cycle_ws = setdiff(all_ws, skip_ws)
# cycle_ws = []

# Go through the loop.
@info "===================== starting ws cycle ========================="
log_error("Diagonalization errors in winding sectors"; init=true)
for ws in cycle_ws
    @info "Winding sector $ws in progress."
    raw_config["winding_sector"] = ws
    try
        multilambda_cycler(param_checks!(raw_config))
    catch LoadError
        log_error("winding_sector: $ws")
        @error "Something went wrong with the diagonalization." ws=ws
    end
end
@info "============================== fin =============================="
