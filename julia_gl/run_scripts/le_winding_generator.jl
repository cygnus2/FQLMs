#===============================================================================

    le_winding_generator.jl - LR, February 2021

    This is a script to create non-zero winding states from maximally flippable
    states in the [0,0,0] winding sector.

===============================================================================#
include("../src/typedefs.jl")
include("../src/io/io.jl")
include("../src/param_checks.jl")
include("../src/io/data_storage.jl")


# Read params and make a logger.
param = param_checks!(read_config("wind_config.yml"); state_run=true)
# param = param_checks!(read_config(ARGS[1]); state_run=true)
latt = LinkLattice(param["L"])

# Set the link type for the correct representation.
const LinkType = latt.S[end]*latt.d >= 63 ? LargeLinkState : SmallLinkState
include("../src/le_stuff/nonzero_wstate_generation.jl")
include("../src/hamiltonian_construction.jl")

function cycle_plaquettes(
        state::LinkType,
        plaquettes::Array{Plaquette,1}
    )::Set{LinkType}
    """ Loops through a list of plaquettes and returns flippable states.
    """
    states = Set{LinkType}()
    for p in plaquettes
        (new_state, _) = apply_u(state, p; sign=false)

        # If U term was not successful, try the U^dagger term.
        # (the order could have been switched - there's always only one
        # possibility for overlap to be generated)
        if isnothing(new_state)
            (new_state, _) = apply_u_dagger(state, p; sign=false)
        end

        # If a state was constructed, we'll add it to the list.
        if !isnothing(new_state)
            push!(states, new_state)
        end
    end
    return states
end


function _exhaust(
        seed::Set{LinkType},
        rest::Set{LinkType},
        level::Integer,
        plaquettes::Array{Plaquette,1},
        max_level::Integer,
        store_states::Bool
    )::Set{LinkType}

    # Some I/O diagnostic stuff.
    @info "Exhausted level $level" num_states=length(seed)
    if store_states
        dump_states(param["state_file"], param["ws_label"]*"/states_lv_$level", collect(seed))
    end

    # ----
    # Start of actual logic.

    # Terminate the recursion after a given # of loops or when there are no
    # new states to be reached.
    if length(seed)==0 || level>=max_level
        @info "Terminating at level $level"
        return union(seed, rest)
    end

    # Loop through all plaquettes and find the flippables.
    # There are several possibilities of doing this, I keep them for the learning
    # effect.

    # 1) Union all sets piece by piece: slow, but memory OK.
    # new_states = Set{LinkType}()
    # for state in seed
    #     new_states = union(new_states, cycle_plaquettes(state, plaquettes))
    # end

    # 2) Union a list of sets. Fast, but limited to a certain length and also
    # potentially memory inefficient.
    # This one is fast, but has a memory problem.
    # l = [cycle_plaquettes(state, plaquettes) for state in seed]
    # new_states = union(l...)

    # 3) Add single pieces to the set - fastest and memory efficient.
    # (actually, BY FAR, the fastest)
    new_states = Set{LinkType}()
    for state in seed
        pstates::Set{LinkType} =  cycle_plaquettes(state, plaquettes)
        for s in pstates
            push!(new_states, s)
        end
    end

    # Prepare for next round and call recursion.
    new_rest = union(seed, rest)
    new_seed = setdiff(new_states, new_rest)
    return _exhaust(new_seed, new_rest, level+1, plaquettes, max_level, store_states)
end


function find_le_states(
        seed_states::Array{LinkType,1},
        latt::LinkLattice;
        max_level::Integer=1000,
        store_states::Bool=true
    )::Array{LinkType,1}

    plaquettes = get_plaquettes(latt)
    new_states = _exhaust(Set{LinkType}(seed_states), Set{LinkType}(), 0, plaquettes, max_level, store_states)
    return collect(new_states)
end


function _parse_base_states(bs::String)::Array{LinkType,1}
    # Remove spaces and split at comma.
    parts = collect(split(filter(x -> !isspace(x), bs[2:end-1]),","))
    return [parse(LinkType, p) for p in parts]
end

# Read the most flippable [0,0,0] states.
new_basis_states = param["base_states"]
if typeof(new_basis_states) == String
    global new_basis_states = _parse_base_states(new_basis_states)
else
    global new_basis_states = LinkType.(param["base_states"])
end

# Loop over directions.
wir = param["winding_sector"]
for (k, dir) in enumerate([xpos, ypos, zpos])
    if wir[k] != 0
        wi = Int(wir[k] / abs(wir[k]))
        for j = wi:wi:wir[k]
            global new_basis_states = increase_winding(new_basis_states, latt, dir; increment=wi)
        end
    end
end

time = @elapsed all_states = find_le_states(new_basis_states, latt; max_level=param["maximum_excitation_level"], store_states=param["store_states"])
@info "Finished finding low energy states." n_states=length(all_states) time=time
