#===============================================================================

    le_winding_generator.jl - LR, February 2021

    This is a script to create non-zero winding states from maximally flippable
    states in the [0,0,0] winding sector. The idea is to only alter these states
    slightly to introduce as litle deviation from the maximally flippable status
    as possible. One way to do this (I'm sure there are others) is to only
    change flip the spins along a straight line between the two sides of the box
    when the line is fully occupied or empty (constant occupancy/spin). For
    periodic boundary conditions, this is essentially a line/loop defect and the
    Gauss law is still intact. It turns out that this produces maximally
    flippable states in the winding sector with one flux along the direction
    of the line.

===============================================================================#
include("../src/typedefs.jl")
include("../src/io/io.jl")
include("../src/io/logging.jl")
include("../src/param_checks.jl")

# Read params and make a logger.
param = read_config("wind_config.yml")
latt = LinkLattice(param["L"])

# Set the link type for the correct representation.
const LinkType = latt.S[end]*latt.d >= 63 ? LargeLinkState : SmallLinkState


function find_lines(latt::LinkLattice, dir::IType)::Array{Array{LinkType,1},1}
    """ Finds the straight lines that connect two sides of the lattice in a given
        direction. These are represented by a list of link indicies.

        For instance on a 2x2x2 lattice we may follow the y-direction starting
        from the lower left corner. The corresponding links would be [2, 8]. Note
        that those are 1-based (Julia default).
    """
    # TODO: implement this.
    return [[5;11]]
end


function increase_winding(state::LinkType, latt::Lattice, dir::IType; increment::IType=1)::Array{LinkType,1}
    """ Takes in a state and returns states with winding changed by one unit
        (positive or negative, toggled by increment) along the specified
        direction (x=1, y=2, z=3).
    """
    lines = find_lines(latt, dir)

    winding_states = Array{LinkType,1}()
    for line in lines
        new_state = state
        for link in line
            new_state = increment > 0 ? create(new_state,link) : annihilate(new_state,link)
            if isnothing(new_state)
                break
            end
        end
        if !isnothing(new_state)
            push!(winding_states, new_state)
        end
    end
    return winding_states
end


# Loop through all provided states. These need to be converted first.
# TODO: what happens for larger dataype - does it work natively? (the reading of the YAML file)
winding_states = Array{LinkType,1}()
base_states = LinkType.(param["base_states"])
for state in base_states
    global winding_states = vcat(winding_states, increase_winding(state, latt, 2; increment=-1))
end

# Check output.
for state in winding_states
    println(Int(state))
end
